#!/usr/bin/env perl
# Copyright (C) 2016-2018 Computomics GmbH

use File::Copy;
use POSIX qw(strftime);
use List::Util qw(min max);
use Config::Simple;

use strict;

my $usage = "$0 <samplefile> <DMR_file> <output_folder> <config file> <context CG|CHG|CHH|combined>\n";

my $samplefile = shift or die $usage;
my $dmrfile = shift or die $usage;
my $output_folder = shift or die $usage;
my $config_file = shift or die $usage;
my $context = shift or die $usage;

die "\n{collectDMRs} ERROR: Cannot find input samplefile\n" if (! -e $samplefile);
die "\n{collectDMRs} ERROR: Cannot find input DMR file\n" if (! -e $dmrfile);
die "\n{collectDMRs} ERROR: Output folder does not exist\n" if (! -e $output_folder);
die "\n{collectDMRs} ERROR: Config file '$config_file' not found!" if (! -e $config_file);

our $CFG;
Config::Simple->new($config_file)->import_names('CFG');


# P-value correction:
print STDERR localtime()." {collectDMRs} p-value correction across $context context\n";
my $ret = system("$CFG::PYTHON_PATH $CFG::SCRIPT_PATH/pv2qv.py $dmrfile $dmrfile.qv");
die "\n{collectDMRs} ERROR: p-value correction failed!\n" if ($ret);

# summarize individual tests into DMR calls and print some useful meta data:
open F, "<$dmrfile.qv" or die "\n{collectDMRs} ERROR: Cannot open $dmrfile.qv\n";


# read sample file for sample order:
my @sample_list;
open S,"<$samplefile" or die "\n{collectDMRs} ERROR: Cannot open $samplefile\n";
while (<S>) {
  my @s = split" ";
  push @sample_list, $s[0];
}
close S;


my $outfile = "$output_folder/DMRs.$context.bed";
if (-e $outfile) {
    move($outfile, "$outfile." . (strftime "%F.%R", localtime));
}

open OUT, ">$outfile" or
  die "\n{collectDMRs} ERROR: Cannot open $outfile for writing\n";

my $lastcoord="";
my %clusters=();
my %context_means=();
my $sites=0;
my $dmr=0; my $hdmr=0;

while (1) {
    my $line = <F>;
    chomp $line;
    my @a = split"\t", $line;

    my $coord = (scalar(@a) == 0) ?
      "" : $a[0]."\t".$a[1]."\t".($a[1]+$a[2]-1)."\t".$a[2]."\t".$a[3];

    if ( $coord ne $lastcoord ) {

        if ($dmr) {
            print OUT $lastcoord . "\t";
            foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                print OUT $cluster.":".$context_means{$cluster} . "\t";
            }
            foreach my $cluster (sort{$a<=>$b} keys %clusters) {
                print OUT $cluster.":".$clusters{$cluster} . "\t";
            }
            print OUT "#:" . $sites . "\t";
            print OUT $context . "\t";
            print OUT $hdmr == 0 ? "-" : $context;
            print OUT "\n";
        }

        %clusters=();
        %context_means=();
        $sites = 0;
        $dmr = 0; $hdmr = 0;

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

    $context_means{$a[4]} = sprintf("%.0f", $a[6]);
    $context_means{$a[5]} = sprintf("%.0f", $a[7]);

    $sites = $a[60];

    $dmr  = 1 if ($a[-1] <= $CFG::FDR_CUTOFF);

    my %min_meth = (
        "CG"       => $CFG::CLUSTER_MIN_METH_CG,
        "CHG"      => $CFG::CLUSTER_MIN_METH_CHG,
        "CHH"      => $CFG::CLUSTER_MIN_METH_CHH,
        "combined" => $CFG::CLUSTER_MIN_METH
    );

    $hdmr = 1 if ($a[-1] <= $CFG::FDR_CUTOFF &&
                  $sites >= $CFG::DMR_MIN_C &&
                  max($a[6], $a[7]) >= $min_meth{$context}/100 &&
                  min($a[6], $a[7]) * $CFG::HDMR_FOLD_CHANGE <= max($a[6], $a[7]));

    $lastcoord = $coord;
}
close F;
close OUT;

system("rm $dmrfile*") if ($CFG::REMOVE_INTMED_FILES);
