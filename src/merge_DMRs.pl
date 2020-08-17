#!/usr/bin/env perl
# Copyright (C) 2016-2018 Computomics GmbH

use File::Copy;
use POSIX qw(strftime);
use List::Util qw(min max);

use strict;

my $usage = "$0 <samplefile> <DMR_file> <output_folder> <python path> <pv2qv.py path> <FDR cutoff> <min. meth. rate> <min. sites> <fold change> <remove intermediate files?>\n";

my $samplefile = shift or die $usage;
my $dmr_file = shift or die $usage;
my $output_folder = shift or die $usage;
my $python_path = shift or die $usage;
my $qv_exec = shift or die $usage;
my $fdr = shift or die $usage;
my $min_meth = shift or die $usage;
my $min_sites = shift or die $usage;
my $fold_change = shift or die $usage;
my $rm_files = shift;

die "\n{mergeDMRs} ERROR: Cannot find input samplefile\n" if (! -e $samplefile);
die "\n{mergeDMRs} ERROR: Cannot find input DMR file\n" if (! -e $dmr_file);
die "\n{mergeDMRs} ERROR: Output folder does not exist\n" if (! -e $output_folder);

die "\n{mergeDMRs} ERROR: Cannot find $qv_exec\n" if (! -e $qv_exec);

# P-value correction:
print STDERR localtime()." {mergeDMRs} p-value correction across all contexts\n";
my $ret = system("$python_path $qv_exec $dmr_file $dmr_file.qv");
die "\n{mergeDMRs} ERROR: p-value correction failed!\n" if ($ret);

print STDERR localtime()." {mergeDMRs} p-value correction for CG context\n";
#system("awk -v OFS='\t' '{print \$9, \$11, \$16, \$10}' $dmr_file > $dmr_file.cg");
$ret = system("
            awk -v OFS='\t' '{print \$9, \$11, \$16, \$10}' $dmr_file > $dmr_file.cg && \
            $python_path $qv_exec $dmr_file.cg $dmr_file.cg.qv");
die "\n{mergeDMRs} ERROR: p-value correction CG failed!\n" if ($ret);

print STDERR localtime()." {mergeDMRs} p-value correction for CHG context\n";
#$awk_cmd = "awk -v OFS='\t' '{print \$26, \$28, \$33, \$27}'";
#$ret = system("$python_path $qv_exec \<($awk_cmd $dmr_file) $dmr_file.qv.chg");
$ret = system("
            awk -v OFS='\t' '{print \$26, \$28, \$33, \$27}' $dmr_file > $dmr_file.chg && \
            $python_path $qv_exec $dmr_file.chg $dmr_file.chg.qv");
die "\n{mergeDMRs} ERROR: p-value correction CHG failed!\n" if ($ret);

print STDERR localtime()." {mergeDMRs} p-value correction for CHH context\n";
#$awk_cmd = "awk -v OFS='\t' '{print \$43, \$45, \$50, \$44}'";
#$ret = system("$python_path $qv_exec \<($awk_cmd $dmr_file) $dmr_file.qv.chh");
$ret = system("
            awk -v OFS='\t' '{print \$43, \$45, \$50, \$44}' $dmr_file > $dmr_file.chh && \
            $python_path $qv_exec $dmr_file.chh $dmr_file.chh.qv");
die "\n{mergeDMRs} ERROR: p-value correction CHH failed!\n" if ($ret);

unlink("$dmr_file") if ($rm_files);
print STDERR localtime()." {mergeDMRs} done\n";


# read sample file for sample order:
my @sample_list;
open S,"<$samplefile" or die "\n{mergeDMRs} ERROR: Cannot open $samplefile\n";
while (<S>) {
  my @s = split" ";
  push @sample_list, $s[0];
}
close S;


# summarize individual tests into DMR calls and print some useful meta data:
prettify_output($dmr_file, $output_folder);


sub prettify_output {
    my ($dmrfile, $output_folder) = @_;

    open F,     "<$dmrfile.qv" or die "\n{mergeDMRs} ERROR: Cannot open $dmrfile.qv\n";
    open F_CG,  "<$dmrfile.cg.qv" or die "\n{mergeDMRs} ERROR: Cannot open $dmrfile.cg.qv\n";
    open F_CHG, "<$dmrfile.chg.qv" or die "\n{mergeDMRs} ERROR: Cannot open $dmrfile.chg.qv\n";
    open F_CHH, "<$dmrfile.chh.qv" or die "\n{mergeDMRs} ERROR: Cannot open $dmrfile.chh.qv\n";

    if (-e "$output_folder/DMRs.bed") {
        move("$output_folder/DMRs.bed", "$output_folder/DMRs.bed.". (strftime "%F.%R", localtime));
    }

    open OUT, ">$output_folder/DMRs.bed" or
      die "\n{mergeDMRs} ERROR: Cannot open $output_folder/DMRs.bed\n";

    open OUT_ALL, ">$output_folder/all_context_DMRs.bed";

    my $lastcoord="";
    my %clusters=();
    my %context_means=();
    my $sites="";
    my %dmr=(); my %hdmr=();
    my $dmr_all=0; my $hdmr_all=0;

    while (1) {
        my $line     = <F>;     chomp $line;     my @a     = split"\t", $line;
        my $line_cg  = <F_CG>;  chomp $line_cg;  my @a_cg  = split"\t", $line_cg;
        my $line_chg = <F_CHG>; chomp $line_chg; my @a_chg = split"\t", $line_chg;
        my $line_chh = <F_CHH>; chomp $line_chh; my @a_chh = split"\t", $line_chh;

        my $coord = (scalar(@a) == 0) ?
          "" : $a[0]."\t".$a[1]."\t".($a[1]+$a[2]-1)."\t".$a[2]."\t".$a[3];

        if ( $coord ne $lastcoord ) {

            if (scalar(keys %dmr) > 0) {
                print OUT $lastcoord . "\t";
                    foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                        print OUT $cluster.":".$context_means{$cluster} . "\t";
                    }
                    foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                        print OUT $cluster.":".$clusters{$cluster} . "\t";
                    }
                    print OUT "#:" . $sites . "\t";
                    print OUT join(",", sort keys %dmr) . "\t";
                    print OUT scalar(keys %hdmr) == 0 ? "-" : join(",", sort keys %hdmr);
                    print OUT "\n";
            }
            if ($dmr_all) {
                print OUT_ALL $lastcoord . "\t";
                foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                    print OUT_ALL $cluster.":".$context_means{$cluster} . "\t";
                }
                foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                    print OUT_ALL $cluster.":".$clusters{$cluster} . "\t";
                }
                print OUT_ALL "#:" . $sites . "\t";
                print OUT_ALL $hdmr_all? "1" : "-";
                print OUT_ALL "\n";
            }

            %clusters=();
            %context_means=();
            $sites = "";
            %dmr = (); %hdmr = ();
            $dmr_all = 0; $hdmr_all = 0;

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

        $context_means{$a[4]} = sprintf("%.0f", $a[6]) . "," .
                                ($a_cg[1] == -1 ? "-" : sprintf("%.0f", $a_cg[1]*100)) . "," .
                                ($a_chg[1] == -1 ? "-" : sprintf("%.0f", $a_chg[1]*100)) . "," .
                                ($a_chh[1] == -1 ? "-" : sprintf("%.0f", $a_chh[1]*100));
        $context_means{$a[5]} = sprintf("%.0f", $a[7]) . "," .
                                ($a_cg[2] == -1 ? "-" : sprintf("%.0f", $a_cg[2]*100)) . "," .
                                ($a_chg[2] == -1 ? "-" : sprintf("%.0f", $a_chg[2]*100)) . "," .
                                ($a_chh[2] == -1 ? "-" : sprintf("%.0f", $a_chh[2]*100));

        $sites = $a[60] . "," . $a_cg[0] . "," . $a_chg[0] . "," . $a_chh[0];

        $dmr_all = 1 if ($a[-1] <= $fdr);
        $dmr{"CG"} = 1 if ($a_cg[-1] <= $fdr);
        $dmr{"CHG"} = 1 if ($a_chg[-1] <= $fdr);
        $dmr{"CHH"} = 1 if ($a_chh[-1] <= $fdr);

        $hdmr_all   = 1 if ($a[-1] <= $fdr &&
                            $a[60] >= $min_sites &&
                            max($a[6], $a[7]) >= $min_meth/100 &&
                            min($a[6], $a[7]) * $fold_change <= max($a[6], $a[7]));
        $hdmr{"CG"} = 1 if ($a_cg[-1] <= $fdr &&
                            $a_cg[0] >= $min_sites &&
                            max($a_cg[1], $a_cg[2]) >= $min_meth/100 &&
                            min($a_cg[1], $a_cg[2]) * $fold_change <= max($a_cg[1], $a_cg[2]));
        $hdmr{"CHG"} = 1 if ($a_chg[-1] <= $fdr &&
                             $a_chg[0] >= $min_sites &&
                             max($a_chg[1], $a_chg[2]) >= $min_meth/100 &&
                             min($a_chg[1], $a_chg[2]) * $fold_change <= max($a_chg[1], $a_chg[2]));
        $hdmr{"CHH"} = 1 if ($a_chh[-1] <= $fdr &&
                             $a_chh[0] >= $min_sites &&
                             max($a_chh[1], $a_chh[2]) >= $min_meth/100 &&
                             min($a_chh[1], $a_chh[2]) * $fold_change <= max($a_chh[1], $a_chh[2]));

        $lastcoord = $coord;
    }
    close F;
    close OUT;
}

system("rm $dmr_file*") if ($rm_files);
