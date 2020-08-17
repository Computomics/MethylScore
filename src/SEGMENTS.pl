#!/usr/bin/env perl
# Copyright (C) 2016-2017 Computomics GmbH

#TODO check sort, python, tabix, python modules

use Getopt::Long;
use Algorithm::Cluster;
use FindBin qw($Bin $RealBin);
use File::Which;
use Data::Dumper;
use List::Util qw(sum max);
use Compress::BGZF::Reader;
use Thread::Pool;
use Thread::Conveyor::Array;
use Thread::Conveyor::Throttled;
use Thread::Conveyor::Tied;
use Thread::Tie::Array;
use File::Copy;
use POSIX qw(strftime);
use Fcntl qw/ :flock /;

use strict;
#use warnings;

my $samplefile; our $matrixfile; my $mrfile;
my $output_folder;
my $OPT_min_covered_sites = 10;
my $OPT_min_coverage = 3;
my $OPT_breakpoint_lookbacks = 0; # NOT USED RIGHT NOW (otherwise set as user parameter)
my $OPT_lookback_maxdistance = 30;
my $OPT_lookback_distance = 30;
my $OPT_MRfreq_fraction = 14;
my $OPT_MR_min_nr_Cs = 5;
my $OPT_window_size;
my $OPT_window_step;
my $OPT_min_cluster_diff = 0;
my $OPT_min_cluster_mean = 0;
my $OPT_fdr = 0.05;
my $OPT_max_lowcov_regionlen = 100;

my $OPT_nr_threads = 1;

my $post_process = 1;

my $OPT_weight_CG = 1;
my $OPT_weight_CHG = 1;
my $OPT_weight_CHH = 1;

my $tabix_exec="";
my $bgzip_exec="";
my $test_binary = "betabin_model";
my $test_exec="$Bin/$test_binary";
my $python_exec="python";
my $qv_exec="";

my $DO_TESTING = 1;
#TODO do something if verbose==1 :-)
my $verbose = 1;
my $DEBUG = 0;


### Get user arguments:
    my %CMD=();
    my $command = GetCom();

    our $OPT_sliding_window = (defined $OPT_window_size && $OPT_window_size > 0) ? 1 : 0;
##################################


### Initialize LOG file:
    my $logfile = "$output_folder/log";
    open LOG, ">>", $logfile;
    print LOG "##########################################\n";
    print LOG "###         Starting DMR calling       ###\n";
    print LOG "##########################################\n";
    print LOG "\n".localtime()."\n";
    print LOG "Executing:\n".$command."\n\n";

    print STDERR "[DMRcalling] ".localtime()." - Starting\n" if ($verbose);

    open CLUSTERTIME, ">$output_folder/clustertime.log";
##################################


### Global variables:
    my @segments=();
    # $ STRING chr
    # $ INT start
    # $ INT end
    # % INT covs  { sampleID => coverage }
    # % INT mrate { sampleID => methrate }
    # $ STRING sample_string
    # $ INT K
    # @ INT cluster_means
    my $Pool = Thread::Pool->new({
        do => 'process_segment',     # this is the function to execute for every job
        workers => $OPT_nr_threads,
        maxjobs => undef,       # setting to undef means: infinite nr of jobs
        #minjobs => $OPT_nr_threads, # job throttling, only spawn jobs if there are this min. nr
        optimize => 'cpu',      # cpu or memory
    });
    my %Jobs=();
    my $Jobs_done=0;
    our %samples=();
    our @sample_list=();
    our %MRpos=();
    #TODO change all sample orders to the variable @sample_list
##################################


### read in MRs of all samples:
    open S,"<$samplefile" or abort("\n{DMRcalling} Cannot open $samplefile\n");
    my $nr_samples=0;
    while (<S>) {
        chomp;
        my @s = split" ";

        my @sample_cols = split ",", $s[1];
        $samples{$s[0]}{"reps"} = \@sample_cols;
        $nr_samples++;
        push @sample_list, $s[0];
        #$samples{$s[0]}{"idx"} = $nr_samples;
        #$samples_order{$s[0]} = $nr_samples;

      	if (-e $s[2] && $mrfile eq "") {
              open F, "<$s[2]" or abort("\n{DMRcalling} ERROR: Cannot open $s[2]\n");
              while (<F>) {

                  chomp;
                  my @a = split" ";

                  if ($a[3] >= $OPT_MR_min_nr_Cs) {
                      push @{$MRpos{$a[0]}{$a[1]}}, \@a;
                      push @{$MRpos{$a[0]}{$a[2]}}, \@a;
                  }

              }
              close F;
      	}

    }
    close S;
    print LOG "".localtime()." - Finished reading sample " . ($mrfile eq "" ? "and MRs " : "") . "file\n";
#print "nr: ".scalar(keys %MRpos)." nrChr1: ".scalar(keys %{$MRpos{"Chr1"}})."\n";
##################################


### read from MRfile already merged across samples, if provided
    if ($mrfile ne "") {
      open F, "<$mrfile" or abort("\n{DMRcalling} ERROR: Cannot open $mrfile\n");
      while (<F>) {
          chomp;
          my @a = split" ";

          if ($a[3] >= $OPT_MR_min_nr_Cs) {
            push @{$MRpos{$a[0]}{$a[1]}}, \@a;
            push @{$MRpos{$a[0]}{$a[2]}}, \@a;
          }
        }
      close F;
      print LOG "".localtime()." - Finished reading MR file\n";
    }
##################################


### Index matrix file if needed:
    if ($DO_TESTING) {

        if ($matrixfile !~ m/.gz$/) {

          if (-e "$matrixfile.processing") {
            warn "{DMRcalling} Matrix file is being processed. Waiting until finished\n";
            print LOG localtime()."Matrix file is being processed. Waiting until finished\n";
          }
          my $max_wait = 1000;
          while (-e "$matrixfile.processing" && $max_wait > 0) { sleep 1; --$max_wait; }
          abort("\n{DMRcalling} ERROR: Compressed matrix not ready.\n") if ($max_wait <= 0);

          if (! -e "$matrixfile.gz" || -z $matrixfile) {
            print LOG localtime()." - Compressing genome matrix file...\n";

            #TODO use BGFZ::Writer!
            system("touch $matrixfile.processing");
            #my $cmd = "$bgzip_exec -f -c $matrixfile > $matrixfile.gz";
            my $cmd = "$bgzip_exec -f $matrixfile";
        	  my $return = system($cmd);
            unlink("$matrixfile.processing");
            abort("\n{DMRcalling} ERROR: Compressing matrix failed\n") if ($return != 0);
            print LOG localtime()." - done\n";
          }

          $matrixfile .= ".gz";
        }

        if (! -e "$matrixfile.tbi") {
            print LOG localtime()." - Indexing genome matrix file...\n";

            #TODO there's also a perl tabix API!
            my $cmd = "$tabix_exec -s 1 -b 2 -e 2 $matrixfile";
            my $return = system($cmd);
            abort("\n{DMRcalling} ERROR: Indexing matrix failed\n") if ($return != 0);
            print LOG localtime()." - done\n";
        }
    }
##################################


### Prepare matrix file, jump to first non-header line:
    my $M = Compress::BGZF::Reader->new_filehandle($matrixfile) or abort("\n{DMRcalling} ERROR: Cannot open matrix file $matrixfile!\n");
    our %M=();
    my $mline = "#";
    while ($mline =~ m/^#/) {
        $mline = <$M>; chomp $mline;
    }
    our @m = split " ", $mline;
    #$M{$m[1]}{"last"} = 0;
##################################


### Output files:
    # don't overwrite previous files:
    if (-e "$output_folder/segments.bed") {
        move("$output_folder/segments.bed", "$output_folder/segments.bed.". (strftime "%F.%R", localtime));
    }
    if (-e "$output_folder/segments.dif") {
        move("$output_folder/segments.dif", "$output_folder/segments.dif.". (strftime "%F.%R", localtime));
    }
##################################




#############################################
######### START ITERATION OVER MR BREAKPOINTS
#############################################
#my @chrs = sort keys %MRpos;

# Get chromosome order:
my $chrs   = `$tabix_exec -l $matrixfile`;
my @chroms = split"\n", $chrs;
foreach my $chr (@chroms) {

    next if (! exists $MRpos{$chr});

    my %MRfreqs=( 0 => 0, );
    my $MRfreq = 0;
    my @positions=();
    my $start_of_MRblock=0;
    my $end_of_MRblock=0;

    foreach my $pos (sort{$a<=>$b} keys %{$MRpos{$chr}}) {
#$DEBUG = 1 if ($pos > 43255577);

        $start_of_MRblock = $pos if ($MRfreq == 0);

        foreach my $MRline (@{$MRpos{$chr}{$pos}}) {

            my @MR = @{$MRline};

            if ($pos == $MR[1]) {
                $MRfreq++;
            }
            elsif ($pos == $MR[2]) {
                $MRfreq--;
            }

        }

        $end_of_MRblock = $pos if ($MRfreq == 0);

        $MRfreqs{$pos} = $MRfreq;
        push @positions, $pos;

print "pos $pos (freq $MRfreq)\n" if ($DEBUG);

        # we now have the MRfreq up to $pos
        # so, we can calculate MR freq differences to previous positions

        # start of methylated block -> start a new segment
        if ($pos == $start_of_MRblock) { #($MRfreqs{$positions[-2]} == 0) {
            my %h = ( "chr" => $chr, "start" => $pos, );
            push @segments, \%h;
print "start of MR block $pos\n" if ($DEBUG);
            next;
        }

        # within or end of MR block:
        else {

            # check if MR frequency across samples changes more than $OPT_MRfreq_fraction
            #   - compare MR freq to previous position, i.e. prev_pos_to_compare = $positions[-2]
            #   - and: compare MR freq to position $OPT_lookback_distance bases upstream of current pos
            #   - or: compare MR freq to $OPT_breakpoint_lookbacks positions before current pos
            #   - if we are at methylated block end: don't compare and take current segment
            #     (prev_pos_to_compare = -1)

            my @prev_pos_to_compare;
            if ($MRfreq == 0) { # this means methylated block ends here!
              push @prev_pos_to_compare, -1;
            }
            else {
print "positions: " . join(",", @positions) . "\n" if ($DEBUG);
              push @prev_pos_to_compare, $positions[-2];

              my $cmp_pos;
              if ($OPT_lookback_distance > 0) {
                my $i=scalar(@positions)-1;
                while ($i>0 && $positions[$i] > $pos-$OPT_lookback_distance) { --$i; }
                $cmp_pos = $positions[$i];
              }
              elsif ($OPT_breakpoint_lookbacks > 0 && scalar(@positions) >= $OPT_breakpoint_lookbacks) {
                # get breakpoint $cmp_pos positions back:
                $cmp_pos = $positions[(-1)*($OPT_breakpoint_lookbacks+1)];
              }

              # limit $cmp_pos to the end of the last segment/start of actual segment:
              $cmp_pos = $segments[-1]{"start"} if (scalar(@segments) > 0 && $segments[-1]{"start"} > $cmp_pos);
              push @prev_pos_to_compare, $cmp_pos if ($cmp_pos != $prev_pos_to_compare[-1]);
            }

print "prev_pos_to_compare: " . join(",", @prev_pos_to_compare) . "\n" if ($DEBUG);

            foreach my $cmp_pos (@prev_pos_to_compare) {

print "pos $pos / freq $MRfreq compared to $cmp_pos / ".$MRfreqs{$cmp_pos}."\n"  if ($DEBUG);

                # if end of meth. block or MR freq change large enough:
                if ( $cmp_pos == -1 ||
                     abs($MRfreq - $MRfreqs{$cmp_pos}) >= $nr_samples * $OPT_MRfreq_fraction/100 ) {

                    my ($adj_start, @subsegments) = get_subsegments($segments[-1]{"chr"}, $segments[-1]{"start"}, $pos);
print "no. subsegments: " . scalar(@subsegments) . "\n" if ($DEBUG);

                    if (scalar(@subsegments) > 0) {
                      ### SUCCESS - there are subsegments with sufficient nr covered sites
                      ### -------

                      pop @segments;
                      foreach my $segm (@subsegments) {
print "######## Final segment: \n" if ($DEBUG);# . %$segm{"chr"}.":".%$segm{"start"}."-".%$segm{"end"}."\n" if ($DEBUG);
                        push @segments, $segm;
                      }

                      # start new (open) segment here unless we're at the MR block end
                      if ($cmp_pos != -1) {
                        my %h = ( "chr" => $chr, "start" => $pos, );
                        push @segments, \%h;
print "OPENING SEGMENT $chr:$pos-\n" if ($DEBUG);
                      }
                    }
                    else {
                      ### NO SUCCESS - upstream segment had insufficient coverage
                      ### ----------

                      # adjust start position if user-defined max. low-cov region length exceeded:
                      $segments[-1]{"start"} = $adj_start;
print "new segment start: $adj_start\n" if ($DEBUG);

                      if ($cmp_pos == -1) {
                        # we're at the end of a MR block and current segment too short -> merge with previous segment - or forget whole MR block:
                        pop @segments;
                        if (scalar(@segments)>0 && $adj_start != $pos) {
                            # overwrite previous segment end and sample coverages:
                            $segments[-1]{"end"} = $pos;
                            (my $p, $segments[-1]{"covs"}, $segments[-1]{"mrate"}) = get_sample_coverages($chr, $segments[-1]{"start"}, $segments[-1]{"end"});
print "UPDATED SEGMENT AT MB end: ".$chr . "\t" . ($segments[-1]{"start"}) . "\t" . ($segments[-1]{"end"}) . "\n" if ($DEBUG);
                        }
                      }
                    }
                } # if MRfreq change large enough

            } # foreach $cmp_pos

        } # within MR block

        # next $pos
        # unless it's the end of an MR block, then:

        # at the end of an MR block, process all segments found within:
        if ($pos == $end_of_MRblock) {

            # cluster segment and do the tests for differential
            # methylation for each pair of clusters:
            foreach my $segm (@segments) {
                if ($OPT_nr_threads > 1) {
                    my $jobid = $Pool->job($segm, \@sample_list, \%samples, $matrixfile);
                    $Jobs{$jobid} = 1;
                }
                else {
                    process_segment($segm, \@sample_list, \%samples, $matrixfile);
                }
            }

            undef @segments;
            undef @positions;
            %MRfreqs=( 0 => 0, );
            %M=();

# if ($pos > 10000000) {
#   print LOG localtime()." Exit at Chr1:10M\n";
#   $Pool->shutdown;
#   exit 0;
# }

        } # end of MR block

#last if ($pos > 10000000);

    } # for each MR breakpoint of a chromosome

    %M=();
    delete $MRpos{$chr};

    print LOG "".localtime()." - Iterated through $chr\n";
    print STDERR "[DMRcalling] ".localtime()." - iterated through $chr\n" if ($verbose);

#last; # skip all Chrs other than Chr1
} # for each chromosome



print LOG localtime()." - Iteration done. Still ".$Pool->workers." threads to finish and ".$Pool->todo." threads pending. Waiting for them...\n";

while (scalar(keys %Jobs) > 0) {
  foreach my $jid (keys %Jobs) {
    my $s = $Pool->result_dontwait($jid);
    if ($s == -1) {
      abort("ERROR: Job $jid failed\n");
    }
    elsif (defined $s) {
      delete $Jobs{$jid};
    }
  }
}


# there might still be some unfinished Jobs
# wait for them:
# $Pool->result further up seems not to wait for the threads to really finish, so:
$Pool->shutdown;


# sort by genomic coords and multiple testing correction:
print LOG localtime()." - All clustering and testing jobs done\n";

if ($post_process) {
  print LOG localtime()." - Now sorting and correcting for multiple testing...\n";
  if (! -e "$output_folder/segments.dif" || -z "$output_folder/segments.dif") {
    die "\n{DMRcalling} ERROR: $output_folder/segments.dif does not exist or is empty!\n";
  }
  my $ret = system("sort -k1,1 -k2,2g -k3,3g $output_folder/segments.dif > $output_folder/segments.dif.sorted");

  die "\n{DMRcalling} ERROR: $output_folder/segments.dif.sorted does not exist or is empty!\n" if ($ret);
  unlink("$output_folder/segments.dif") if (!$ret);

  $ret = system("$python_exec $qv_exec $output_folder/segments.dif.sorted $output_folder/segments.dif.sorted.qv");
  die "\n{DMRcalling} ERROR: $output_folder/segments.dif.sorted.qv does not exist or is empty!\n" if ($ret);
  unlink("$output_folder/segments.dif.sorted") if (!$ret);
  print LOG localtime()." - done\n" if (!$ret);

  # summarize individual tests into DMR calls and print some useful meta data:
  prettify_output("$output_folder/segments.dif.sorted.qv");
  print LOG localtime()." - DMRs written into $output_folder/DMRs.bed\n";

  print LOG "\n".localtime()." - DMRcalling FINISHED!\n";
}
else {
  print LOG "\n".localtime()." - FINISHED\n";
  print STDERR "[DMRcalling] ".localtime()." - FINISHED\n";
}
close LOG;


exit 0;




#################################
############ Subroutines
#################################


# This routine is run in parallel for each segment:
sub process_segment {
    my ($segment_unclustered, $smplist, $smpls, $matrixfile) = @_;
    my @sample_list = @{$smplist};
    my %samples = %{$smpls};

    # cluster samples on the particular segment:
    # (cluster function also checks and ignores samples with too few covered sites)
    my $segment_clustered = cluster_segment($segment_unclustered);
    my %segm = %{$segment_clustered};
    undef $segment_unclustered;
    undef $segment_clustered;

#print $segm{"chr"} . "\t" . $segm{"start"} . "\t" . $segm{"end"} . "\t" . $segm{"K"} . "\n";

    # write segment into file independent of whether samples have been grouped (K>1) or not (K=1):
print "PRINT TO SEGMENTS.BED: ".$segm{"start"}."-".$segm{"end"}."\n" if ($DEBUG);
    open OUTPUT_SEGMENTS, ">>$output_folder/segments.bed";
    flock(OUTPUT_SEGMENTS, LOCK_EX);
    print OUTPUT_SEGMENTS   $segm{"chr"} . "\t" .
                            $segm{"start"} . "\t" .
                            ($segm{"end"}-$segm{"start"}) . "\t" .
                            $segm{"sample_string"} . "\n";
    close OUTPUT_SEGMENTS;


#TODO: extend to SGE submission system, or offer Ramdisk solution!
#TODO include flag if one wants to test for variability even if there are no clusters found! (i.e. if K=1)

    # TEST the clusters against each other:
    if ($segm{"K"} > 1 && max(@{$segm{"cluster_means"}}) >= $OPT_min_cluster_mean && $DO_TESTING) {

        my $nr = int(rand(100000));
        my $filename = "$output_folder/TMPbbtest$nr." . $segm{"chr"} . "-" . $segm{"start"} . "-" . $segm{"end"};

        # for the statistical test we need the matrix file, and to drastically speed up the betabin program we extract only the region of the segment (otherwise the whole huge matrix is read in each instance):
        my $return = system(
            "$tabix_exec $matrixfile " .
            $segm{"chr"}.":".$segm{"start"}."-".$segm{"end"} .
            "> $filename.matrix 2>> $logfile"
        );
        return -1 if ($return != 0);

        if (-z "$filename.matrix") {
            print "\n{DMRcalling} ERROR: Matrix empty! $filename.matrix @ " .
            $segm{"chr"}.":".$segm{"start"}."..".$segm{"end"} . "\n";
            return -1;
        }

        # we need to gather the columns of the samples in $cid1 and $cid2 and their replicates
        # this is needed for the differential test as input
        # columns of samples in the matrix have been stored in %samples
        my @sample_cols = ();
        for (my $i=0; $i<length($segm{"sample_string"}); ++$i) {
            #push @sample_cols, "";
            if (substr($segm{"sample_string"}, $i, 1) =~ /(\d)/) {
                my $cid = $1-1;
                my $sample_name = $sample_list[$i];
                foreach my $col (@{$samples{$sample_name}{"reps"}}) {
                    $sample_cols[$cid] .= "," if (defined $sample_cols[$cid] && $sample_cols[$cid] ne "");
                    $sample_cols[$cid] .= $col-4;
                }
#print "sample cols of cid $cid: ".$sample_cols[$cid]."\n";
            }
        }

        # do pairwise comparisons between clusters:
        my $file_input;
        for (my $cid1 = 0; $cid1<$segm{"K"}-1; $cid1++) {
            for (my $cid2=$cid1+1; $cid2<$segm{"K"}; $cid2++) {
#print "[$coords] cmps: $cid1 - $cid2\n";
                # create input file of betabin_model to test for differential methylation:
                $file_input = $filename . "-$cid1-$cid2";
                open TMP, ">$file_input";
                print TMP   $segm{"chr"} . "\t" .                       # chr
                            $segm{"start"} . "\t" .                     # start pos
                            ($segm{"end"} - $segm{"start"}) . "\t" .    # end pos
                            $segm{"sample_string"} . "\t" .             # segment length
                            ($cid1+1) . "\t" .                          # cluster1 ID
                            ($cid2+1) . "\t" .                          # cluster2 ID
                            $segm{"cluster_means"}[$cid1] . "\t" .      # cluster1 mean meth.
                            $segm{"cluster_means"}[$cid2] . "\t" .      # cluster2 mean meth.
                            "\n";
                close TMP;

                # bash script to execute betabin_model instance:
                my $file_exec = $file_input . ".sh";
                open TMP, ">$file_exec"
                    or print "\n{DMRcalling} ERROR: Cannot create $file_exec\n" and return -1;
                print TMP "$test_exec difftest $filename.matrix $file_input " .
                          $sample_cols[$cid1] . "\\\|" . $sample_cols[$cid2] .
                          " $file_input.dif 1 0.001 6 &> /dev/null\n";
                close TMP;

#print "now: $file_exec\n" if ($DEBUG);

                # Finally, do the test for differential methylation between cid1 and cid2:
#print "[$coords] Run test...\n";
                my $cmd = "bash $file_exec";
                system("$cmd");
#print "[$coords] done\n";

                # read out the result and write it to the file ALL_TESTS:
                my $file_prefix = $file_exec;
                $file_prefix =~ s/.sh//;

                if (! -e "$file_prefix.dif") {
                    print "\n{DMRcalling} ERROR: BetaBin Test failed for $file_exec\n";
                    return -1;
                }
                if (-z "$file_prefix.dif") {
                    print "\n{DMRcalling} ERROR: dif file empty! Unprocessed: $file_exec. Continue anyway.\n";
                    return -1;
                }

#print "[$coords] open $file_prefix.dif and $output_folder/segments.dif\n";
                open TMP, "<$file_prefix.dif"
                    or print "\n{DMRcalling} ERROR: Cannot open $file_prefix.dif\n" and return -1;
                open ALL_TESTS, ">>$output_folder/segments.dif"
                    or print "\n{DMRcalling} ERROR: Cannot open $output_folder/segments.dif\n" and return -1;
                flock(ALL_TESTS, LOCK_EX)
                    or print "\n{DMRcalling} ERROR: Cannot lock segments.dif\n" and return -1;
                print ALL_TESTS <TMP>
                    or print "\n{DMRcalling} ERROR: Couldn't write test result to $file_prefix.dif.\n" and return -1;
                close ALL_TESTS;
                close TMP;

                # remove all files for that specific cluster comparison:
                unlink($file_input)
                    or print "\n{DMRcalling} WARNING: Cannot remove $file_input\n";
                unlink($file_exec)
                    or print "\n{DMRcalling} WARNING: Cannot remove $file_exec\n";
                unlink("$file_input.dif")
                    or print "\n{DMRcalling} WARNING: Cannot remove $file_input.dif\n";
                # foreach (glob "$file_prefix*") {
                #     unlink($_) or print LOG "WARNING: Cannot delete $_!\n";
                # }

            }
        } # pairwise cluster comparisons

        # remove matrix file for that specific segment that was processed in this function:
        unlink("$filename.matrix") or print "\n{DMRcalling} WARNING: Cannot remove $filename\n";

    } # segment tested (K>1)

}


sub get_sample_coverages {
    my ($chr, $start, $end) = @_;
#print "REGION $start - $end\n" if ($DEBUG);
    my %PASS = ();
    my %sites_per_sample=();
    my %rates_per_sample=();

    my $pos = $start;

    if (!exists $M{$pos}) {
        # browse up to $pos:
#print "m0: " . $m[0] . " m1: " . $m[1] . " pos $pos\n" if ($DEBUG);
        while ($chr ne $m[0] || $m[1] < $pos) {
            $mline = <$M>; chomp $mline;
            abort("Chromosome $chr not found in matrix!") if ($mline eq "");
            @m=split/\t/, $mline;
        }
    }
print "  GET_COVERAGES, start-end $start-$end pos $pos m1 ".$m[1]."\n" if ($DEBUG);
    if (defined $OPT_window_size && !exists $M{$pos} && $pos < $m[1]) {
        # due to the sliding window analysis, there can be start positions that are not a cytosine,
        # so jump to the next one:
        while (!exists $M{$pos} && $pos < $m[1]) { $pos++; print "m1/pos: $m[1]/$pos\n" if ($DEBUG); }
    }

    while ($pos < $end) {

        #foreach my $sample (sort{$samples_order{$a}<=>$samples_order{$b}} keys %samples_order) {
        foreach my $sample (@sample_list) {

            # initialize $sites_per_sample and $rates_per_sample
            $sites_per_sample{$sample} = 0 if (!exists $sites_per_sample{$sample});
            $rates_per_sample{$sample} = [] if (!exists $rates_per_sample{$sample});

            if (!exists $M{$pos}{$sample}) {
                if ($m[1] != $pos) { die "m1 $m[1] unequal pos $pos, sample $sample, start $start end $end, segment-1: ".$segments[-1]{"start"}." segments-2: ".$segments[-2]{"start"}."-".$segments[-2]{"end"}.", segments-3: ".$segments[-3]{"start"}."-".$segments[-3]{"end"}."\n"; }

                # read the actual matrix content:
                my @read_counts=(0,0);
                foreach my $rep_col (@{$samples{$sample}{"reps"}}) {
                    my @tmp = (0,0,0,0);
                    @tmp = split"/", $m[$rep_col] if ($m[$rep_col] ne "0");
                    $read_counts[0] += $tmp[2];
                    $read_counts[1] += $tmp[3];
                }

#print "  pos $pos: $sample cov ".($read_counts[0]+$read_counts[1])."\n" if ($DEBUG);

                my $cov = $read_counts[0]+$read_counts[1];

                if ($cov >= $OPT_min_coverage) {
                    my $rate = 100*$read_counts[0]/$cov;

                    # weight the rate by context:
                    $rate = ($rate * $OPT_weight_CG > 100) ? 100 : $rate * $OPT_weight_CG if ($m[2] eq "CG");
                    $rate = ($rate * $OPT_weight_CHH > 100) ? 100 : $rate * $OPT_weight_CHH if ($m[2] eq "CHH");
                    $rate = ($rate * $OPT_weight_CHG > 100) ? 100 : $rate * $OPT_weight_CHG if ($m[2] eq "CHG");

                    $M{$pos}{$sample}{"cov"} = 1;
                    $M{$pos}{$sample}{"mrate"} = $rate;
                    $sites_per_sample{$sample}++;
                    push @{$rates_per_sample{$sample}}, $rate;
                }
                else {
                    $M{$pos}{$sample}{"cov"} = 0;
                    $M{$pos}{$sample}{"mrate"} = -1;
                    push @{$rates_per_sample{$sample}}, -1;
                }
#print "  sample $sample mrate ".$M{$pos}{$sample}{"mrate"}."\n" if ($DEBUG);
            }
            else {
                $sites_per_sample{$sample} += $M{$pos}{$sample}{"cov"};
                push @{$rates_per_sample{$sample}}, $M{$pos}{$sample}{"mrate"};
            }
#print "cov sample $sample: " . $M{$pos}{$sample}{"cov"} . "\n" if ($DEBUG);

            $PASS{$sample} = 1 if ($sites_per_sample{$sample} >= $OPT_min_covered_sites);
        }

        # iterate further:
        if (!exists $M{$pos}{"next"}) {
            $mline = <$M>; chomp $mline;
            @m=split/\t/, $mline;
            last if ($mline eq "" || $m[0] ne $chr);
            $M{$pos}{"next"} = $m[1];
            #$M{$m[1]}{"last"} = $pos;
            $pos = $m[1];
#print "  new m1:  ".$m[1]."\n" if ($DEBUG);
        }
        else {
            $pos = $M{$pos}{"next"};
#print "  new pos: $pos\n" if ($DEBUG);
        }
    }
print "   -> pass ".(%PASS && scalar(keys %PASS)>1 ? 1 : 0)."\n" if ($DEBUG);

    return (
        (%PASS && scalar(keys %PASS) > 1 ? 1 : 0),
        \%sites_per_sample,
        \%rates_per_sample
    );
}


sub get_subsegments {
    my ($chr, $start, $end) = @_;

print "SEGMENT $start-$end\n" if ($DEBUG);

    if (! $OPT_sliding_window) {
      # make sure that while-loop further down is executed
      # only once if sliding window is turned off:
      $OPT_window_size = $end - $start;
      $OPT_window_step = $end - $start;
    }

    my @subsegment = ();
    my $adj_start = $start;
    my $global_end = $end;
    $end = $start + $OPT_window_size;
    $end = $global_end if ($end > $global_end); # limit to end of segment, of course

    while (1) {#($end < $global_end) {
        my ($pass, $sites_per_sample, $rates_per_sample);
        ($pass, $sites_per_sample, $rates_per_sample) = get_sample_coverages($chr, $start, $end);
print "  subsegment covs: " . Dumper($sites_per_sample) if ($DEBUG);

        if ($pass) {
            # subsegment passes coverage filter -> store it
            my %h = ( "chr" => $chr, "start" => $start, "end" => $end, "covs" => $sites_per_sample, "mrate" => $rates_per_sample );
            push @subsegment, \%h;

            $start += $OPT_window_step;
        }
        else {
            # get second highest no of covered sites across samples
            # if that's below MIN_COV_SITES_PER_SAMPLE, don't store subsegment
            my $max2nd_cov_sites_per_samples = (sort values %$sites_per_sample)[-2];

            if ($max2nd_cov_sites_per_samples == 0 ||
                $end - $start > $OPT_max_lowcov_regionlen) {
              # try new segment from position $end on:
              $start = $end;
              $end = ($end + $OPT_window_size) > $global_end ? $global_end : ($end + $OPT_window_size);
              $adj_start = $start;
print "new start: $start\n" if ($DEBUG);
              next unless ($start == $end); # directly go to get_sample_coverages at the beginning of while-loop
            }
            elsif (scalar(@subsegment) > 0) {
                # update previous subsegment if current subsegment has insufficient coverage, but more than $CFG::MIN_COV_SITES_PER_SAMPLE:
                $subsegment[-1]{"end"} = $end;
                ($pass, $subsegment[-1]{"covs"}, $subsegment[-1]{"mrate"}) =
                  get_sample_coverages($chr, $subsegment[-1]{"start"}, $end);
                $start += $OPT_window_step;
            }
        }
        last if ($end == $global_end);
        $end = ($end + $OPT_window_step > $global_end) ? $global_end : $end + $OPT_window_step;
    }

    # discard subsegments that have not been closed (e.g. if only 1 sample had coverage at end of segment)
    #pop @subsegment if (! defined $subsegment[-1]{"end"});

    return ($adj_start, @subsegment);
}


# deprecated
sub sliding_window {
    my ($s) = @_;
    my %segment = %{$s};

print "segment ".$segment{"start"}."-".$segment{"end"}."\n" if ($DEBUG);

    my @subsegment=();

    my $start = $segment{"start"};
    my $global_end = $segment{"end"};
    my $end = $start + $OPT_window_size;
    while ($end < $global_end) {
        my ($pass, $sites_per_sample, $rates_per_sample);
        #($pass, $segments[-1]{"covs"}, $segments[-1]{"mrate"}) = get_sample_coverages($chr, $segments[-1]{"start"}, $pos);
        ($pass, $sites_per_sample, $rates_per_sample) = get_sample_coverages($segment{"chr"}, $start, $end);
        if ($pass) {
            my %h = ( "chr" => $segment{"chr"}, "start" => $start, "end" => $end, "covs" => %{$sites_per_sample}, "mrate" => %{$rates_per_sample} );
            push @subsegment, \%h;
            $start += $OPT_window_step;
        }
        elsif (scalar(@subsegment) > 0) {
            $subsegment[-1]{"end"} = $end;
            ($pass, $subsegment[-1]{"covs"}, $subsegment[-1]{"mrate"}) = get_sample_coverages($segment{"chr"}, $subsegment[-1]{"start"}, $end);
            $start += $OPT_window_step;
        }
        $end = ($end + $OPT_window_step > $global_end) ? $global_end : $end + $OPT_window_step;
    }

    return \@subsegment;
}


sub cluster_segment {
    my ($s) = @_;
    my %segment = %{$s};

my $start_time = time();
my $mid_time = undef;
print "cluster ".$segment{"start"}."-".$segment{"end"}."\n" if ($DEBUG);

    my (@matrix, @mask, @samples_covered);
    my (@cluster_means, @cluster_errors);
    my @sample_means;
    my $sample_idx = 0;
    my %sample_string=();

    my ($rate_cumsum , $rate_occ) = 0;
    foreach my $sample (sort keys %{$segment{"mrate"}}) {

        # disregard samples that have less covered sites than $OPT_min_covered_sites
        if ($segment{"covs"}{$sample} < $OPT_min_covered_sites) {
            if ($segment{"covs"}{$sample} == 0) { $sample_string{1} .= "."; }
            else                                { $sample_string{1} .= "-"; }
            next;
        }
        else                                    { $sample_string{1} .= "o"; }

        #my @tmp=();
        @sample_means = [];#\@tmp;

        # fill @matrix:
        my @mrate_vec = @{$segment{"mrate"}{$sample}};
#print "  sample $sample cov ".$segment{"covs"}{$sample}." rates ".join(",", @mrate_vec)."\n" if ($DEBUG);
        push @matrix, \@mrate_vec;

        # fill @mask:
        my @mask_vec = @mrate_vec;
        for (my $i=0; $i<scalar(@mrate_vec); ++$i) {
            if ($mrate_vec[$i] < 0) {
                $mask_vec[$i] = 0;
            }
            else {
                $mask_vec[$i] = 1;

                $rate_cumsum += $mrate_vec[$i];
                $rate_occ++;

                push @{$sample_means[$sample_idx]}, $mrate_vec[$i];
            }
        }
        push @mask, \@mask_vec;

        # print sample mean
print "--> sample mean: ".(sum(@{$sample_means[$sample_idx]})/scalar(@{$sample_means[$sample_idx]}))."\n" if ($DEBUG);
        ++$sample_idx;
    }

    my @globalmean = [ $rate_cumsum / $rate_occ ];
    push @cluster_means, \@globalmean;
print "==> cluster mean: ".($rate_cumsum / $rate_occ)."\n" if ($DEBUG);

$mid_time = time() - $start_time;

    # Do the CLUSTERING with different K (1<K<5):
    my $K=2; my @Wks=(); my @sse=(); my %clusters=();
    for (; $K <= int(100/$OPT_min_cluster_diff); ++$K) {

        # break if there are less covered samples than $K:
        last if (scalar(@matrix) < $K);

        my %params = (
            nclusters => $K,
            data => \@matrix,
            mask => \@mask,
            transpose => 0,
            npass => 100,
            method => 'a',
            dist => 'e',
            initialid => [],
        );
        #TODO maybe test to use max and min sample as initial cluster centers? (maybe only for half of the iterations and check with the other random centers)
        #TODO parameterize maxK and maxIterations of clustering

        my ($clusters, $error, $found) = Algorithm::Cluster::kcluster(%params);

        ####################################################
        # Calculate cluster means and evaluate their differences to each other
        ####################################################

print "K=$K - error $error - found $found\n" if ($DEBUG);
print Dumper($clusters)."\n" if ($DEBUG);

        # first calculate cluster means and store all data points for each cluster (cluster_points):
        $sample_idx = 0; my %rate_sums=(); my %cluster_points=();
        foreach my $sample (sort keys %{$segment{"mrate"}}) {

            if ($segment{"covs"}{$sample} < $OPT_min_covered_sites) {
                if ($segment{"covs"}{$sample} == 0) { $sample_string{$K} .= "."; }
                else                                { $sample_string{$K} .= "-"; }
                next;
            }

            my $clusterID = @{$clusters}[$sample_idx];

            $sample_string{$K} .= ($clusterID+1);

            if (!exists $cluster_points{$clusterID}) {
                my @tmp=(); $cluster_points{$clusterID} = \@tmp;
            }

            foreach my $rate (@{$segment{"mrate"}{$sample}}) {
                if ($rate >= 0) {
                    push @{$cluster_points{$clusterID}}, $rate;
                    $rate_sums{$clusterID} += $rate; #sum( @{$segment{"mrate"}{$sample}} );
                }
            }

            $sample_idx++;

        }

        # evaluate differences between cluster means and discard K if differences too low:
        my @means=(); my $PASS=1;
        for (my $clusterID=0; $clusterID<$K; ++$clusterID) {
            push @means, ( $rate_sums{$clusterID} / scalar(@{$cluster_points{$clusterID}}) );
print "mean cluster $clusterID: ".$means[$clusterID]."\n" if ($DEBUG);
            for (my $c=0; $c<$clusterID; ++$c) {
print "cluster diffs $c $clusterID: ".abs($means[$c] - $means[$clusterID])."\n" if ($DEBUG);
                $PASS = 0 if ( abs($means[$c] - $means[$clusterID]) < $OPT_min_cluster_diff );
            }
        }

        if (! $PASS) {
            last;
        }

        push @cluster_means, \@means;

    } # for each K


#        # calculate compactness of clustering acc. to "Elbow" heuristic
#
#         # then calculate Dk's and Wk by iterating over all cluster points
#         # (notation according to
#         # https://datasciencelab.wordpress.com/2013/12/27/finding-the-k-in-k-means-clustering/)
#         my %means=(); my $Wk = 0; my $sse=0;
#         for (my $clusterID=0; $clusterID<$K; ++$clusterID) {
#             # foreach CLUSTER
#
#             # calculate cluster mean:
#             $means{$clusterID} = $rate_sums{$clusterID} / scalar(@{$cluster_points{$clusterID}});
#
#             #my $Dk=0;
#             foreach my $rate (@{$cluster_points{$clusterID}}) {
#                 # foreach METHRATE
#
#                 #$Dk += abs($rate - $means{$clusterID});
#                 $sse += ( ( $rate - $means{$clusterID} ) * ( $rate - $means{$clusterID} ) );
#             }
#
#             #$Wk += sqrt( $Dk );
#         }
#
# print "SSE diff to prev K: ".($sse[-1] - $sse)."\n" if ($DEBUG);
#         push @sse, $sse;
#
#         # my $kk=2;
#         # foreach my $wk (@Wks) {
#         #     print "Wk diff to K=$kk: ".($Wk-$wk)."\n" if ($DEBUG);
#         #     $kk++;
#         # }
#         # push @Wks, $Wk;
#
# #print "Wk $Wk\n" if ($DEBUG);
# foreach (keys %means) { print "mean cluster $_: ".$means{$_}."\n" if ($DEBUG); }
# print "--------\n" if ($DEBUG);
#
#     }

    --$K;
print "#### K = $K\n" if ($DEBUG);
my $total_time = time() - $start_time;
print CLUSTERTIME $mid_time . "\t" . $total_time . "\n";

    $segment{"K"} = $K;
    $segment{"sample_string"} = $sample_string{$K};
    $segment{"cluster_means"} = $cluster_means[-1];

    return \%segment;
}


sub prettify_output {
    my ($dmr_file) = @_;

    open F, "<$dmr_file" or die "\n{DMRcalling} ERROR: Cannot open $dmr_file\n";
    if (-e "$output_folder/DMRs.bed") { move("$output_folder/DMRs.bed", "$output_folder/DMRs.bed.". (strftime "%F.%R", localtime)); }
    open OUT, ">$output_folder/DMRs.bed";

    my $lastcoord="";
    my %clusters=();
    my %clustermeans=();
    my $dmr=0;

    while (1) {
        my $line = <F>; chomp $line;
        my @a = split"\t", $line;

        my $coord = (scalar(@a) == 0) ? "" : join("\t", @a[0..3]);

        if ( $coord ne $lastcoord ) {

            if ($dmr) {
                print OUT $lastcoord . "\t";
                    foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                        print OUT $cluster.":".$clustermeans{$cluster} . "\t";
                    }
                    foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                        print OUT $cluster.":".$clusters{$cluster} . "\t";
                    }
                    print OUT "\n";
            }

            %clusters=();
            %clustermeans=();
            $dmr = 0;

        }

        last if ($coord eq "");

        if (scalar(keys %clusters) == 0) {
            my @clusterstring = split"", $a[3];

            my $i=0;
            foreach my $clusterID (@clusterstring) {
                if ($clusterID ne "-" && $clusterID ne ".") {
                    $clusters{$clusterID} .= "," if (exists $clusters{$clusterID});
                    $clusters{$clusterID} .= $sample_list[$i];
                }
                $i++;
            }
        }

        $clustermeans{$a[4]} = sprintf("%.0f", $a[6]);
        $clustermeans{$a[5]} = sprintf("%.0f", $a[7]);

        $dmr = 1 if ($a[-1] <= $OPT_fdr);

        $lastcoord = $coord;
    }
    close F;
    close OUT;
}


sub abort {
  my ($str) = @_;

  if ($OPT_nr_threads > 1) {
    $Pool->shutdown;
  }

  die $str;
}


### Read command line parameters --------------------------------------------------------------
sub GetCom {

	my @usage = ("$0

Select segments to test for differential methylation based on methylated regions (MRs)
of a population of samples

Mandatory parameters:
-s   FILE           Sample file (3-column file containing:
                      * sampleID
                      * comma-separated list of columns in genome matrix file
                        containing sample's replicates, 0-based (i.e. specify
                        '4' for sample 1 in column 5 of matrix)
                      * sample's MR file)
-m   FILE           Genome Matrix file
-o   STRING         Output folder

Optional parameters:
          default:
-v   INT        3   Minimum coverage at each cytosine sites within segment to test
-n   INT       10   Minimum number of cytosine sites with min.coverage >= c within
                     segment to test (in at least two samples)
-a   INT      100   Maximal length (bp) of low-coverage region (not fulfilling -v
                     and -n criteria); segment will be discarded if exceeding
-p   DOUBLE    14   Percentage of samples that have to change methylation status at
                     a given position to delimit segment to test
-d   INT       30   Position <d> bp upstream to which MR frequency change is compared
-i   INT        0   Minimum cluster mean difference between any pair of clusters to
                     perform statistical test
-j   INT       20   Minimum mean methylation level of at least one cluster
-f   DOUBLE  0.05   Maximum false discovery rate (FDR)

-w   INT      off   Sliding window size for long segments
-x   INT      off   Sliding window step size (maximum is <w>)

-z   INT        1   Number of threads

-B   STRING         Path of the statistical test executable
                      \"$test_binary\". Default: $Bin/$test_binary
-T   STRING         Path of the 'tabix' executable. Default: 'which tabix'
-E   STRING         Path of the 'bgzip' executable. Default: 'which bgzip'
-Y   STRING         Path of the pv2qv.py script for multiple testing correction.
                     Default: $Bin
-K   STRING         Python executable (default: 'python')

-q                  Quiet mode

Advanced parameters (experimental):
-l                  Skip testing for differential methylation, only retrieve segments
-r   FILE           Already merged MR file
-c   DOUBLE     1   Weight for methylation level in CG context
-g   DOUBLE     1   Weight for methylation level in CHG context
-h   DOUBLE     1   Weight for methylation level in CHH context

Copyright (C) 2016-2017, Computomics GmbH
\n");

	die(@usage) if (@ARGV == 0);

	my $command = $0 . " " . join(" ", @ARGV);

	GetOptions(\%CMD, "s=s", "m=s", "o=s", "c=s", "v=i", "n=i", "a=i", "p=f", "d=i", "w=i", "x=i", "z=i", "cg=f", "chg=f", "chh=f", "i=f", "j=i", "f=f", "B=s", "T=s", "E=s", "Y=s", "K=s", "l", "r=s", "q", "u=i", "no-post-process");

	die("Please specify strain file!\n[ABORTED]\n") unless defined($CMD{s});
	die("Please specify matrix file!\n[ABORTED]\n") unless defined($CMD{m});
	die("Please specify output folder!\n[ABORTED]\n") unless defined($CMD{o});

	$samplefile = $CMD{s};
	$matrixfile = $CMD{m};
	$output_folder = $CMD{o};
	chop $output_folder if ($output_folder =~ m"\/$");
  mkdir $output_folder if (! -e $output_folder);

  if (defined($CMD{v})) { $OPT_min_coverage = $CMD{v}; }
  if (defined($CMD{n})) { $OPT_min_covered_sites = $CMD{n}; }
  if (defined($CMD{a})) { $OPT_max_lowcov_regionlen = $CMD{a}; }
  if (defined($CMD{p})) { $OPT_MRfreq_fraction = $CMD{p}; }
  if ($OPT_MRfreq_fraction > 0 && $OPT_MRfreq_fraction < 1) {
    die "{DMRcalling} Parameter MR frequency change -p must be an integer (0<p<100)\n";
  }
  if (defined($CMD{d})) { $OPT_lookback_distance = $CMD{d}; }
  if (defined($CMD{w})) { $OPT_window_size = $CMD{w}; }
  if (defined($CMD{x})) { $OPT_window_step = $CMD{x}; }
  if (defined($CMD{"cg"})) { $OPT_weight_CG = $CMD{"cg"}; }
  if (defined($CMD{"chg"})) { $OPT_weight_CHG = $CMD{"chg"}; }
  if (defined($CMD{"chh"})) { $OPT_weight_CHH = $CMD{"chh"}; }
  if (defined($CMD{i})) { $OPT_min_cluster_diff = $CMD{i}; }
  if ($OPT_min_cluster_diff > 0 && $OPT_min_cluster_diff < 1) {
    die "{DMRcalling} Parameter min. cluster difference -i must be an integer (0<i<100)\n";
  }
  if (defined($CMD{j})) { $OPT_min_cluster_mean = $CMD{j}; }
  if ($OPT_min_cluster_mean < 0 || $OPT_min_cluster_mean > 100) {
    die "{DMRcalling} Parameter min. cluster mean must be between 0 and 100\n";
  }
  if (defined($CMD{f})) { $OPT_fdr = $CMD{f}; }
  if ($OPT_fdr < 0 || $OPT_fdr > 1) {
    die "{DMRcalling} Parameter FDR must be between 0 and 1\n";
  }
  if ($OPT_fdr > 0.5) { warn "{DMRcalling} Good idea to set FDR as high as $OPT_fdr?\n"; }
  if (defined($CMD{z})) { $OPT_nr_threads = $CMD{z}; }

  if (defined($CMD{"no-post-process"})) { $post_process = 0; }

  if (defined($CMD{r})) { $mrfile = $CMD{r}; }

  if (defined($CMD{u})) { $DEBUG = 1 if ($CMD{u}==1); }

  $OPT_window_step = $OPT_window_size if (defined $OPT_window_size &&
                                       (! defined $OPT_window_step || $OPT_window_step == 0));
  if (defined $OPT_window_size && defined $OPT_window_step && $OPT_window_step > $OPT_window_size) {
    $OPT_window_step = $OPT_window_size;
    warn "{DMRcalling} Sliding window step not allowed to be larger than window size, set to window size\n";
  }

  ###### CHECKING TABIX AND BGZIP:
  if (defined($CMD{T})) { $tabix_exec = $CMD{T}; }
	if ($tabix_exec eq "") {
		$tabix_exec = which('tabix');
		if ($tabix_exec eq "") {
        		die "\n{DMRcalling} ERROR!!! Program 'tabix' not found! Please install or specify its path via -T!\n[ABORTED]\n";
		}
	}
  if (defined($CMD{E})) { $bgzip_exec = $CMD{E}; }
      if ($bgzip_exec eq "") {
              $bgzip_exec = which('bgzip');
              if ($bgzip_exec eq "") {
                      die "\n{DMRcalling} ERROR!!! Program 'bgzip' not found! Please install or specify its path via -E!\n[ABORTED]\n";
              }
      }
  if (defined($CMD{K})) { $python_exec = $CMD{K}; }


  ###### CHECKING BETABIN_MODEL:
  if (defined($CMD{B})) { $test_exec = $CMD{B}; }
  if (! -e "$test_exec") {
      die "\n{DMRcalling} Cannot find statistical test binary '$test_binary'!\n";
  }

  ###### CHECKING pv2qv:
  if (defined($CMD{Y})) { $qv_exec = $CMD{Y}; }
  elsif (-e "$Bin/pv2qv.py") { $qv_exec = "$Bin/pv2qv.py"; }
  else { die "\n{DMRcalling} Cannot find Storey's test binary pv2qv.py (looked in $Bin)!\n"; }

  if (defined($CMD{l})) { $DO_TESTING = 0; }
  if (defined($CMD{q})) { $verbose = 0; }

  return $command;
	#print LOG "\n#".localtime()."\n# Selecting regions to test\n\nExecuting:\n".$command."\n";
}
