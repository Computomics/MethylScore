#!/usr/bin/env perl
#  Copyright (C) 2013-14 Joerg Hagmann, Max Planck Institute for 
#				Developmental Biology, Tuebingen
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

use Getopt::Long;

my $samplefile;
my $reffile;
my $iformat="mx";
my $outfile;
my $mergedfile="";
my $not_core = 0;
my $header = 1;
my $tmpdir = ".";
my $verbose = 1;

my %CMD;
GetCom();

# read reference for determining methylation class later on:
my %REF=();
my $chr="";
my $seq="";
if ($iformat !~ m/^mx/) {
open R,"<$reffile" or die "Cannot open $reffile\n";
while( <R> ) {
	chomp;
	if($_ =~ /^>/ ) {
		if($seq ne "") {
			$REF{$chr} = $seq;
		}
		$chr = substr($_,1);
                # deal with fasta headers containing more than one word
                $chr = (split(' ',$chr))[0];
		$seq = "";
	}
	else {
		$seq .= $_;
	}
}
# make sure bases are uppercase
$REF{$chr} = uc $seq;
close R;

print STDERR "reference has been read\n" if ($verbose);
}

# get all strain data and give strains in a fixed order (%strains)
my %strains=();
my @files=();
my $c=0;


# store header
my $header_str="";
$header_str .= "#chr\tpos\tclass\tstrand" if ($header);

#BRR57	/ebio/abt6_read_storage/backup/data/solexa_analysis/ATH/Genome/Haplotype1/AlignmentFolder/BRR57/MethylAnalysis02/ConsensusAnalysis/methylated_sites.txt
open F,"<$samplefile" or die "Cannot open $samplefile\n";
while (<F>) {
	chomp;
	next if ($_ =~ m/^#/);

	my @s = split/\t/,$_;

	if (! -e $s[1]) { die "Cannot find $s[1]!\n"; }

	if (($iformat eq "bismark" || $iformat eq "mx") && $mergedfile eq "") {
                # print sample name into bismark input file:
	        open FILE, "<$s[1]";
       	        my $rand = int(rand(100000));
	        open OUT, ">$s[1].tmp$rand";
        	while (my $line = <FILE>) {
                	print OUT $s[0]."\t".$line;
               	}
               	close FILE;
               	close OUT;

	        push @files, "$s[1].tmp$rand";
        }
        elsif ($iformat eq "shore" || $iformat eq "mxX") {
		push @files, $s[1];
	}

	$strains{$s[0]} = $c;
	$c++;

	# extend header
	$header_str .= "\t$s[0]" if ($header);
}
close F;

my $nr_strains = $c;


# merge meth sites of specified strains
my $nr = int(rand(100000));
if ($mergedfile eq "") {
	print STDERR "merging and sorting methylation consensus files...\n";

	while (-e "tmp$nr.all_methsites.txt") {
		$nr = int(rand(100000));
	}

	my $cmd = "sort -m -k2,2d -k3,3g -T $tmpdir ";
	foreach (@files) {
		$cmd .= "$_";
		$cmd .= " " if ($_ ne $files[scalar(@files)-1]);
	}
	$cmd .= " > $tmpdir/tmp$nr.all_methsites.txt";

	print STDERR $cmd."\n";
	system("$cmd");
	$mergedfile = "$tmpdir/tmp$nr.all_methsites.txt";
}
#}
#else {
#	print STDERR "Skip merging methylated_sites files - tmp_all_methsites_$str.txt already exists\n";
#}


my @cols = (8,12,14);

#two different formats:
#0	1	2	3	4	5	6	7	8	9	10	11	12	13	14
#met1-3	1	110	G	9	0/9	0	9	0/0	1	4/5	4/5	0/5	0	40	

#0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	
#6G0086  Chr1    1095    C       1       0/1     0/0     0/0     0       0       0/0     REP     1       BAL     1/0     1/0     CL      0/1     0/0     1       1       1       1       0       14      core_support_hard, pc_hard,

# parse through merged file and create matrix
print STDERR "parsing merged file..\n" if ($verbose);

open OUT, ">$outfile" or die "Cannot open $outfile!\n";
print OUT $header_str."\n" if ($header);

my $skipped_sites = 0;

open A,"<$mergedfile" or die "Cannot open $mergedfile\n";
$chr="";
$pos=0;
my @strand=();
my %context=();
my @a=();
my %line=();
my $discard = 0;
while (1) {
	my $line = <A>; chomp $line;
	next if ($line =~ m/^#/);
	@a=split/\t/, $line;

	if (scalar(@a) > 15) {
		if ($not_core) {
			@cols = (10,17,24);
		}
		else {
			@cols = (10,18,24);
		}
	}

	if ($pos != 0 && ($a[1] ne $chr || $a[2] > $pos+2)) {
		# if there's sth to print:
		foreach my $prev_pos (sort{$a<=>$b} keys %line) {

			if ($prev_pos < $pos || $pos eq "") {
				# get methylation class first:
				$discard = 0;
				my $class = "?";
				if ($strand{$prev_pos} eq "C") {
					# plus strand
					my $triplet;
					if ($iformat eq "shore" || $iformat eq "bismark") {
						$triplet = substr($REF{$chr}, $prev_pos-1, 3);
					} elsif ($iformat eq "mx") {
						$triplet = $context{$prev_pos};
					}

					if ($iformat eq "mxX") {
						$class = $context{$prev_pos};
					}
					else {
						if ($triplet =~ m/^CG/) {
							$class = "CG";
						}
						elsif ($triplet =~ m/^C.G/) {
							$class = "CHG";
						}
						elsif ($triplet =~ m/^C/) {
							$class = "CHH";
						}
						elsif ($iformat !~ /^mx/ && length($triplet) < 3) {
							$discard = 1;
						}
						else {
							print STDERR "Warning: Sequence context not recognized: $triplet at $chr:$prev_pos in reference file\n" if ($verbose);
							$skipped_sites++;
							$discard = 1;
						}
					}
				}
				else {
					# minus strand
					my $triplet;
					if ($iformat eq "shore" || $iformat eq "bismark") {
						$triplet = substr($REF{$chr}, $prev_pos-3, 3);
					} elsif ($iformat eq "mx") {
						$triplet = $context{$prev_pos};
					}

					if ($iformat eq "mxX") {
						$class = $context{$prev_pos};
					}
					else {
						if ($triplet =~ m/CG$/) {
							$class = "CG";
						}
						elsif ($triplet =~ m/^C.G/) {
							$class = "CHG";
						}
						elsif ($triplet =~ m/G$/) {
							$class = "CHH";
						}
						elsif ($iformat !~ m/^mx/ && length($triplet) < 3) {
							$discard = 1;
						}
						else {
							print STDERR "Warning: Sequence context not recognized: $triplet at $chr:$prev_pos in reference file\n" if ($verbose);
							$skipped_sites++;
							$discard = 1;
						}
					}
				}
	
				if (scalar(keys %{$line{$prev_pos}}) > 0 && !$discard) {
					print OUT $chr."\t".$prev_pos."\t".$class."\t".$strand{$prev_pos};
	
					# print whole genome matrix row:
					for (my $i=0; $i!=$nr_strains; ++$i) {
						if (exists $line{$prev_pos}{$i}) {
							print OUT "\t".$line{$prev_pos}{$i};
						}
						else {
							print OUT "\t0";
						}
					}
					print OUT "\n";
				}

				delete $line{$prev_pos};
				delete $strand{$prev_pos};
				delete $context{$prev_pos} if ($iformat =~ /^mx/);
			}
		}

		# empty %line for next entry
		#%line = ();
	}

	if ($line eq "") {
                # END OF PROGRAM:
                close A;
		print STDERR "\n$skipped_sites skipped sites due to ambiguous contexts (commonly containing N's),\n";
		print STDERR "\tor because of incomplete contexts (at scaffold/chromosome borders.\n";
		unlink("$tmpdir/tmp$nr.all_methsites.txt");
                #print STDERR "NOTE: Did not rm $tmpdir/tmp$nr.all_methsites.txt\n### FINISHED\n";
                exit;
        }

# mx output:
#6G0086EPI       Chr1    1222    CG      5       5       37      2       13      37
#6G0086EPI       Chr1    1222    CGG     .       .       .       2       15      37
#6G0086EPI       Chr1    1223    GGG     .       .       .       2       13      38

	my @pos=();

	$chr = $a[1];
	$pos = $a[2];
	if ($iformat eq "shore") {
		$strand{$pos} = $a[3];
	}
	elsif ($iformat eq "bismark") {
                $pos = $a[3];
                $strand{$pos} = substr($REF{$chr}, $pos-1, 1);

		if ($strand{$pos} ne "C" && $strand{$pos} ne "G") {
	                # IGNORE genome positions where there is no C or G !!!!!!
			$skipped_sites++;
			print STDERR "Warning: position $chr:$pos does not contain C or G in reference! Skipped\n" if ($verbose);
			next;
	        }
        }

	if ($iformat eq "mx") {
		if ($a[3] =~ m/[^ACGT]/) { # discard heterozygous sites:
			$skipped_sites++; next;
		}

		#if ($a[3] eq "CGG") { # special case CGG
		#	$strand{$pos} = "G";
		#	$pos += 2;
		#}
		#els
		if ($a[3] =~ m/^C/ && $a[4] ne ".") { # triplet starts with C and is covered
			$strand{$pos} = "C";
		}
		elsif ($a[4] eq ".") { # look for next unprocessed G of that sample
			$pos = $pos + index($a[3], "G", $pos-$a[2]) - ($pos-$a[2]);
			while (exists $line{$pos}{$strains{$a[0]}}) {
				$pos = $pos + index($a[3], "G", $pos-$a[2]+1) - ($pos-$a[2]);
			}
			$strand{$pos} = "G";
			#}
			#else { # no unique G in context -> heterozygous site -> discard
			#	$skipped_sites++; next;
			#}
		}
		$context{$pos} = $a[3];
		#else { # all other cases, where there are no unique C's or G's in context string
		#	$skipped_sites++; next;
		#}
	}

	if ($iformat eq "mxX") {
		#if ($a[3] =~ m/[^.ACGT]/) { # discard heterozygous sites:
		#	$skipped_sites++; next;
		#}

		if ($a[4] ne ".") { # triplet starts with C and is covered
			if ($a[3] !~ m/^C/) { $skipped_sites++; next; } # discard heterozygous/SNP sites
			$strand{$pos} = "C";
		}
		elsif ($a[4] eq ".") {
			my $p = ($a[3] =~ m/^CG\./ ? 1 : 2);
			$pos+=$p;
			if (substr($a[3], index($a[3], "\.")+1+$p, 1) ne "G") { $skipped_sites++; next; } # discard het. sites
			$strand{$pos} = "G";
		}
		$context{$pos} = substr($a[3], 0, index($a[3], "\."));
	}


#BRR57	1	22	C	1	1/0	0/0	0/0	0	0	0/0	REP	1	BAL	1/0	1/0	CL	1/0	0/0	1	1	0	0	0	14	core_support_hard, clones_hard,
#0-8-86	1	32	C	1	1/0	1/0	100/100	1	0	100/100	REP	1	BAL	1/0	0/1	CL	1/0	1/0	0	0	1	1	0	14	core_support_hard, clones_hard,

	# fill %line
	if (!$discard) {
		my $qual; my @s; my $rate; my @support;
                if ($iformat eq "shore") {
			$qual = $a[$cols[2]];
			@s = split/\//, $a[$cols[0]]; # 8/10
			@support = split/\//, $a[$cols[1]]; #12/18	# cn cc_support
			if ($not_core) {
				$rate = sprintf("%.1f", $support[0]*100/($support[0]+$support[1]));
			} else {
				$rate = sprintf("%.1f", $s[1]);	# cn cc_meth_rate
			}
		}
		elsif ($iformat eq "bismark") {
                        $qual = 40;
                        $rate = sprintf("%.1f", $a[4]);
                        @support = ($a[5], $a[6]);
                }
                elsif ($iformat =~ /^mx/) {
                        $qual = ($strand{$pos} eq "C" ? $a[6] : $a[9]);
                        $rate = ($strand{$pos} eq "C" ? sprintf("%.1f", $a[4]*100/$a[5]) : sprintf("%.1f", $a[7]*100/$a[8]));
                        @support = ($strand{$pos} eq "C" ? ($a[4], $a[5]-$a[4]) : ($a[7], $a[8]-$a[7]));
                }

		if ($support[0]+$support[1] > 0) {
			$line{$pos}{$strains{$a[0]}} = $qual."/".$rate."/".$support[0]."/".$support[1]; # 14/24
		}
	}

	if ($iformat =~ /^mx/ && $a[4] ne "." && $a[7] ne "." && index($a[3], "G") >= 0) { # additional position?
		#$pos = $pos + index($a[3], "G", $pos-$a[2]+1) - ($pos-$a[2]);
		my $context = ($iformat eq "mx" ? $a[3] : substr($a[3], 0, index($a[3], "\.")));
		my $p = ($context eq "CG" ? 1 : 2);
		$pos += $p;
		$strand{$pos} = "G";
		if (substr($a[3], index($a[3], "\.")+1+$p, 1) ne "G") { $skipped_sites++; next; } # discard het. sites
		$qual = $a[9];
		$rate = sprintf("%.1f", $a[7]*100/$a[8]);
		@support = ($a[7], $a[8]-$a[7]);
		if ($support[0]+$support[1] > 0) {
			$line{$pos}{$strains{$a[0]}} = $qual."/".$rate."/".$support[0]."/".$support[1];
			$context{$pos} = $context;
		}
	}
}
close A;
close OUT;




sub GetCom {

	my @usage = ("$0

Generates raw genome matrix from several samples

-s	STRING	Sample file (2-column file containing:
			* sampleID
			* path of methylation consensus file generated
			  by the program specified with -i
-i	STRING	Input format ('shore', 'bismark', 'mx', 'mxX', default: $iformat)
		mx: methylExtract, mxX: methylExtract with regular seq. contexts
-r	STRING	Genome reference fasta file (not needed for 'mx'/'mxX')
-o	STRING	Output file
-f	STRING	Merged (over all samples) and coordinate-sorted
		methylation consensus file (sorting step will be skipped)
-C		Use entire read sequence, i.e. do not restrict to core region
		(only for 'shore')
-h		Disable header in genome matrix containing sample IDs
-t	STRING	temporary directory (default: .)

-q		Quiet mode

Copyright (C) 2013-15 Joerg Hagmann,
Max Planck Institute for Developmental Biology, Tuebingen, Germany.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License v3 or later.
\n");

	die(@usage) if (@ARGV==0);
	GetOptions(\%CMD, "s=s", "i=s", "o=s", "r=s", "f=s", "C", "h", "t=s", "q");

	die("Please specify sample file!\n[ABORTED]\n") unless defined($CMD{s});
	die("Please specify output file!\n[ABORTED]\n") unless defined($CMD{o});

	$samplefile = $CMD{s};
	$reffile = $CMD{r};
	$outfile = $CMD{o};

	if (defined $CMD{i}) { $iformat = $CMD{i}; }
	if (defined $CMD{o}) { $outfile = $CMD{o}; }
	if (defined $CMD{f}) { $mergedfile = $CMD{f}; }
	if (defined $CMD{C}) { $not_core = 1; }
	if (defined $CMD{h}) { $header = 0; }
	if (defined $CMD{q}) { $verbose = 0; }
	if (defined $CMD{t}) { $tmpdir = $CMD{t}; }

	die("Please specify reference file!\n[ABORTED]\n") unless defined($CMD{r} || $iformat =~ m/^mx/);
}


