#!/usr/bin/env perl
# Copyright (C) 2016-2018 Computomics GmbH

use strict;

my $usage = "$0 <samplefile> <outfile> <nr. MR blocks in file>\n";

my $samplefile = shift or die $usage;
my $outfile = shift or die $usage;
my $batch_MRblock = shift or die $usage;

my %MRfreqs;
my @samples;

my @CHR;
my $CHR_FILLED = 0;

# read in sample data and sample-specific MR files:
open S, "<$samplefile";
while (<S>) {

  chomp;
  my @s = split" ";

  push @samples, $s[0] . "\t" . $s[1];

  if (-e $s[2]) {
				my $prev_chr = "";
				open F, "<$s[2]" or die "\n{split_MRfile} ERROR: Cannot open $s[2]";
				while (<F>) {
				
					chomp;
					my @a = split" ";
					my ($chr, $start, $end) = @a;
					push @a, $s[0];

					push @{$MRfreqs{$chr}{$start}}, \@a;
					push @{$MRfreqs{$chr}{$end}}, \@a;

					if ($chr ne $prev_chr && !$CHR_FILLED) {
							push @CHR, $chr;
					}

					$prev_chr = $chr;
				}
				close F;
				$CHR_FILLED = 1;
	}
}
close S;



# iterate over MR breakpoints and output MRs in chunks
# of <batch_MRblock> MR blocks:
my $nr_MRblocks = 0;
my $MRfreq = 0;

my $offset = 0;

open OUT, ">$outfile.0";

foreach my $chr (@CHR) {
#foreach my $chr (sort keys %MRfreqs) {

	foreach my $pos (sort{$a<=>$b} keys %{$MRfreqs{$chr}}) {

		foreach my $MRline (@{$MRfreqs{$chr}{$pos}}) {

				my @MR = @{$MRline};
		
				if ($pos == $MR[1]) {
						$MRfreq++;
						print OUT join("\t", @MR) . "\n";
				}
				elsif ($pos == $MR[2]) {
						$MRfreq--;
				}
		}

		if ($MRfreq == 0) {
			$nr_MRblocks++;

			if ($nr_MRblocks % $batch_MRblock == 0) {
				close OUT;

				my $batch_nr = int($nr_MRblocks / $batch_MRblock) + $offset;
				open OUT, ">$outfile." . $batch_nr;
			}
		}
	}
	close OUT;

	my $batch_nr = int($nr_MRblocks / $batch_MRblock) + 1;
	$offset++;
	open OUT, ">$outfile." . $batch_nr;
}

close OUT;