#!/usr/bin/perl -w

=about
###########################################################################
###########################################################################

  *****************************************  
  *****  MethylExtract (version 1.9.1) ****  
  *****************************************  


 Computational Epigenomics and Bioinformatics
  Dept. of Genetics & Inst. of Biotechnology 
             University of Granada           

         Web: http://bioinfo2.ugr.es/        

 This program is Copyright (C) 2012-16:
 Ricardo Lebr�n (rlebron@ugr.es), Guillermo Barturen (bartg01@gmail.com), Antonio Rueda (aruemar@gmail.com), Jos� L. Oliver (oliver@ugr.es), Michael Hackenberg (hackenberg@ugr.es)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 For questions, feedback, etc. please contact to:
 Ricardo Lebr�n (rlebron@ugr.es), Guillermo Barturen (bartg01@gmail.com), Jos� L. Oliver (oliver@ugr.es), Michael Hackenberg (mlhack@gmail.com)

 To see the options, please launch MethylExtract without any command line arguments


############################################################################
############################################################################
=cut

use strict;
use threads;
use threads::shared;
use Thread::Queue;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);

#################
#Default options#
#################
my %optionsDefault = (
	seq => "NA", inDir => "NA", outDir => "NA",
	qscore => "phred33-quals", p => 4,
	minQ => 20, methNonCpGs => 0.9, varFraction => 0.1,
	context => "CG", maxPval => 0.05, minDepthMeth => 1, maxStrandBias => 0.7,
	bedOut => "N", wigOut => "N", delDup => "N",
	FirstIgnor => 0, LastIgnor => 0,
	minDepthSNV => 1,
	flagW => "NA", flagC => "NA", 
	tagW => "NA", tagC => "NA", #(deprecated input parameters)
	simDupPb => 32, chromDiv => 400, memNumReads => 200000, 
	chromSplitted => "N", 
	chromSorted => "N", #(deprecated input parameter)
	peOverlap => "N",
	samtools => "samtools"
);
#################
#################

our $bedOut=$optionsDefault{bedOut};
our $wigOut=$optionsDefault{wigOut};
our $vcfOut=$optionsDefault{vcfOut};

my $discardPos: shared=0;
my $discardReads: shared=0;
my $duplicatedReadsDel: shared=0;
my $discardDepthPos: shared=0;
my $uncheckedDepthSNV: shared=0;
my $strandBiasDiscards: shared=0;

my (%homoSNV,%heteroSNV): shared;
my (%numContext): shared;
my @chroms: shared;
my $time_start=time;

###############  Getting options  #############

my ($inDir,$seq,$context,$checkBisulfite,$methNonCpGs,$qscore,$Q,$minSNVperc,$outDir,$p,$delDup,
	$FirstIgnor,$tagW,$tagC,$simDupPb,$chromDiv,$memNumReads,$chromSort,$minDepthSNV,$maxPval,
	$minDepthMeth,$maxStrandBias,$typeFasta,$multiSeqID,$peOverlap,$LastIgnor,$samtools);
my %subQscore=("phred33-quals"=>\&sanger,"solexa-quals"=>\&solexa,"phred64-quals"=>\&illumina);
my %bases=("AG"=>"R","GA"=>"R","CT"=>"Y","TC"=>"Y","GC"=>"S","CG"=>"S","AT"=>"W","TA"=>"W","GT"=>"K"
			,"TG"=>"K","AC"=>"M","CA"=>"M","A"=>"A","C"=>"C","G"=>"G","T"=>"T","CGT"=>"B","CTG"=>"B"
			,"GCT"=>"B","GTC"=>"B","TCG"=>"B","TGC"=>"B","ACG"=>"V","AGC"=>"V","CAG"=>"V","CGA"=>"V"
			,"GCA"=>"V","GAC"=>"V","ACT"=>"H","ATC"=>"H","TCA"=>"H","TAC"=>"H","CTA"=>"H","CAT"=>"H"
			,"AGT"=>"D","ATG"=>"D","GTA"=>"D","GAT"=>"D","TGA"=>"D","TAG"=>"D","N"=>"N");
my %CGseq=("CG"=>"A");
my %CHGseq=("CAG"=>"A","CTG"=>"A","CCG"=>"D","CGG"=>"R");
my %CHHseq=("CAA"=>"D","CAT"=>"D","CAC"=>"D","TTG"=>"R","ATG"=>"R","GTG"=>"R",
			"CTA"=>"D","CTT"=>"D","CTC"=>"D","TAG"=>"R","AAG"=>"R","GAG"=>"R",
			"CCA"=>"D","CCT"=>"D","CCC"=>"D","TGG"=>"R","AGG"=>"R","GGG"=>"R");
			
($inDir,$seq,$context,$checkBisulfite,$methNonCpGs,$qscore,$Q,$minSNVperc,$outDir,$bedOut,
	$wigOut,$p,$delDup,$FirstIgnor,$tagW,$tagC,$vcfOut,$simDupPb,$chromDiv,$memNumReads,$chromSort,
	$minDepthSNV,$maxPval,$minDepthMeth,$maxStrandBias,$typeFasta,$multiSeqID,$peOverlap,
	$LastIgnor,$samtools)=
&GetOptions($context,$checkBisulfite,$methNonCpGs,$qscore,$Q,$minSNVperc,$outDir,$bedOut,$wigOut,
	$p,$delDup,$FirstIgnor,$tagW,$tagC,$vcfOut,$simDupPb,$chromDiv,$memNumReads,$chromSort,
	$minDepthSNV,$maxPval,$minDepthMeth,$maxStrandBias,$typeFasta,$multiSeqID,$peOverlap,
	$LastIgnor,$samtools);

##############   Starting MethylExtract  ############

print "\n################### Running MethylExtract v1.9 ###################\n\n";

print "\nFunctions summary:\n";
if ($typeFasta eq "Single") {print "Single-fasta references in $seq\n"}
else {print "Multi-fasta reference selected in $seq\n"}

my ($tagWout,$tagCout);
foreach my $keys (keys %{$tagW}) {$tagWout.="$keys".","}
chop($tagWout);
foreach my $keys (keys %{$tagC}) {$tagCout.="$keys".","}
chop($tagCout);
print "Watson strand FLAG: $tagWout & Crick strand FLAG: $tagCout\n";
if ($delDup eq "N") {print "Duplicated reads will not be deleted (default)\n"}
else {
	print "Duplicated reads will be deleted, seed differences must be higher than $simDupPb\n";
	print "Warning: The de-duplication step is not advised for RRBS methodology\n";
}
if ($methNonCpGs==0) {print "Deleting bisulfite unconverted reads has been deactivated\n"}
if ($FirstIgnor==0) {print "5' bases of the reads will not be ignored (default)\n"}
if ($LastIgnor==0) {print "3' bases of the reads will not be ignored (default)\n"}
if ($peOverlap eq "N") {print "Pair-end overlapping parts will be included (default)\n"}
else {print "Pair-end mate-2 overlapping parts will be discarded\n"}
if ($maxStrandBias==0) {print "SNVs strand bias threshold has been deactivated\n"}
print "samtools path: $samtools\n";

print "\nQuality values resume:\n";
print "Selected Qscore: $qscore\n";
print "Bisulfite check fraction/integer: $methNonCpGs\n";
print "Bases from the 5' end to be ignored: $FirstIgnor\n";
print "Bases from the 3' end to be ignored: $LastIgnor\n";
print "Minimum PHRED score for methylation and base call: $Q\n";
print "Minimum depth for methylation calls: $minDepthMeth\n";
print "Minimum depth for SNV calls: $minDepthSNV\n";
print "Alleles minimum frequency to be considered: $minSNVperc\n";
print "Maximum allowed strand bias: $maxStrandBias\n";
print "SNV p-value threshold: $maxPval\n\n";

print "Activated Outputs:\n";
print "Default methylation levels output: $context context\n";
if ($bedOut eq "Y") {print "BED output format\n"}
if ($wigOut eq "Y") {print "WIG output format\n"}
if ($minSNVperc==1) {print "SNVs checker has been deactivated\n";}
elsif ($minSNVperc>0.5) {print "For varFraction>0.5 SNVs outputs will not be shown, try lower varFraction to extract SNVs\n"}
else {print "SNVs VCF output format\n";}
print "\n";

#### Split alignments in chromosomes #############

my $QueueSplit=Thread::Queue->new();
my $ResultsSplit=Thread::Queue->new();
if ($chromSort eq "N") {
	print "Sorting alignments by chromosome\n";
	my $returned=&StaticthreadsSplit($QueueSplit,$ResultsSplit,$inDir,\&splitChroms,$seq);
}
else {
	my %test=();
	if ($typeFasta eq "Single") {
		foreach (glob "$seq/*.{fa,fas,fasta,fna}") {
			$_=~/(\S+)\.(fa|fas|fasta|fna)$/i;
			my @splitFile=split(/\//,$1);
			$test{$splitFile[$#splitFile]}=1;
		}
	}
	elsif ($typeFasta eq "Multi") {
		foreach (@{$multiSeqID}) {$test{$_}=1;}
	}
	else {}
	foreach (glob "$inDir/*.{sam,sam.gz,bam}") {
		$_=~/(\S+)\.(sam|sam.gz|bam)$/i;
		my @splitFile=split(/\//,$1);
		if (!$test{$splitFile[$#splitFile]}) {die "Input files don't seem to be grouped by chromosome, please select chromSplitted=N\n"}
	}
}

############## MAIN PROCESS ################

print "Starting main process\n";
my $Queue=Thread::Queue->new();
my $Results=Thread::Queue->new();
my $returned=&Staticthreads($Queue,$Results,$seq,\&Main);

############     OUTPUT FILES    ################

my $coverContextCG_ref=\my %coverContextCG;
my $coverContextCHG_ref=\my %coverContextCHG;
my $coverContextCHH_ref=\my %coverContextCHH;
my $countMeth_ref=\my %countMeth;
my $countUnmeth_ref=\my %countUnmeth;
my $sumTotal_ref=\my %sumTotal;
print "\nMerging Methylation & SNVs Outputs\n";
#Outputs
if ($context eq "CG" or $context eq "ALL") {
	open OUT, ">$outDir/CG.output";
	print OUT "#CHROM\tPOS\tCONTEXT\tWatson METH\tWatson COVERAGE\tWatson QUAL\tCrick METH\tCrick COVERAGE\tCrick QUAL\n";
	for (my $j=0;$j<=$#chroms;$j++) {
		if (-e "$outDir/CG_$chroms[$j].output") {
			open IN, "$outDir/CG_$chroms[$j].output";
			while (my $line=<IN>) {
				print OUT "$line";
				($coverContextCG_ref,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref)=&countValues($line,"CG",$coverContextCG_ref,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
			}
			close IN;
		}
		unlink ("$outDir/CG_$chroms[$j].output");
	}
	close OUT;
}
if ($context eq "CHG" or $context eq "ALL") {
	open OUT, ">$outDir/CHG.output";
	print OUT "#CHROM\tPOS\tCONTEXT\tWatson METH\tWatson COVERAGE\tWatson QUAL\tCrick METH\tCrick COVERAGE\tCrick QUAL\n";
	for (my $j=0;$j<=$#chroms;$j++) {
		if (-e "$outDir/CHG_$chroms[$j].output") {
			open IN, "$outDir/CHG_$chroms[$j].output";
			while (my $line=<IN>) {
				print OUT "$line";
				($coverContextCHG_ref,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref)=&countValues($line,"CHG",$coverContextCHG_ref,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
			}
			close IN;
		}
		unlink ("$outDir/CHG_$chroms[$j].output");
	}
	close OUT;
}
if ($context eq "CHH" or $context eq "ALL") {
	open OUT, ">$outDir/CHH.output";
	print OUT "#CHROM\tPOS\tCONTEXT\tWatson METH\tWatson COVERAGE\tWatson QUAL\tCrick METH\tCrick COVERAGE\tCrick QUAL\n";
	for (my $j=0;$j<=$#chroms;$j++) {
		if (-e "$outDir/CHH_$chroms[$j].output") {
			open IN, "$outDir/CHH_$chroms[$j].output";
			while (my $line=<IN>) {
				print OUT "$line";
				($coverContextCHH_ref,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref)=&countValues($line,"CHH",$coverContextCHH_ref,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
			}
			close IN;
		}
		unlink ("$outDir/CHH_$chroms[$j].output");
	}
	close OUT;
}

if ($minSNVperc<=0.5) {
	open VCF, ">$outDir/SNVs.vcf";
	print VCF "##fileformat=VCFv4.1\n##source=MethylExtract 1.9\n##reference=folder:$seq\n##phasing=unphased\n";
	print VCF "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n";
	print VCF "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Reads depth (unaffected nucleotides by bisulfite conversion)\">\n";
	print VCF "##INFO=<ID=MQ,Number=.,Type=Float,Description=\"Qscore average\">\n";
	print VCF "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias (1-Fisher's exact test p-value)\">\n";
	print VCF "##FORMAT=<ID=BQ,Number=.,Type=Float,Description=\"Average base quality for reads supporting alleles. For each allele, in the same order as listed\">\n";
	print VCF "##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">\n";
	print VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA0001\n";
	for (my $j=0;$j<=$#chroms;$j++) {
		if (-e "$outDir/SNVs_$chroms[$j].output") {
			open IN, "$outDir/SNVs_$chroms[$j].output";
			while (my $line=<IN>) {print VCF "$line"}
			close IN;
		}
		unlink ("$outDir/SNVs_$chroms[$j].output");
	}
	close VCF;
}
else {}

if ($bedOut eq "Y") {
	if ($context eq "CG" or $context eq "ALL") {
		open BEDOUT, ">$outDir/CG.bed";
		print BEDOUT "track name=CG Methylation  description=MethylExtract_1.9 visibility=1 useScore=1\n";
		for (my $j=0;$j<=$#chroms;$j++) {
			open IN, "$outDir/CG_$chroms[$j].bed";
			while (my $line=<IN>) {print BEDOUT "$line"}
			close IN;
			unlink ("$outDir/CG_$chroms[$j].bed");
		}
		close BEDOUT;
	}
	if ($context eq "CHG" or $context eq "ALL") {
		open BEDOUT, ">$outDir/CHG.bed";
		print BEDOUT "track name=CHG Methylation  description=MethylExtract_1.9 visibility=1 useScore=1\n";
		for (my $j=0;$j<=$#chroms;$j++) {
			open IN, "$outDir/CHG_$chroms[$j].bed";
			while (my $line=<IN>) {print BEDOUT "$line"}
			close IN;
			unlink ("$outDir/CHG_$chroms[$j].bed");
		}
		close BEDOUT;
	}
	if ($context eq "CHH" or $context eq "ALL") {
		open BEDOUT, ">$outDir/CHH.bed";
		print BEDOUT "track name=CHH Methylation  description=MethylExtract_1.9 visibility=1 useScore=1\n";
		for (my $j=0;$j<=$#chroms;$j++) {
			open IN, "$outDir/CHH_$chroms[$j].bed";
			while (my $line=<IN>) {print BEDOUT "$line"}
			close IN;
			unlink ("$outDir/CHH_$chroms[$j].bed");
		}
		close BEDOUT;
	}
}

if ($wigOut eq "Y") {
	if ($context eq "CG" or $context eq "ALL") {
		open WIGOUT, ">$outDir/CG.wig";
		print WIGOUT "track type=wiggle_0 name=CG Methylation  description=MethylExtract_1.9 glyph=wiggle_density bgcolor=red fgcolor=red\n";
		for (my $j=0;$j<=$#chroms;$j++) {
			open IN, "$outDir/CG_$chroms[$j].wig";
			while (my $line=<IN>) {print WIGOUT "$line"}
			close IN;
			unlink ("$outDir/CG_$chroms[$j].wig");
		}
		close WIGOUT;
	}
	if ($context eq "CHG" or $context eq "ALL") {
		open WIGOUT, ">$outDir/CHG.wig";
		print WIGOUT "track type=wiggle_0 name=CHG Methylation  description=MethylExtract_1.9 glyph=wiggle_density bgcolor=red fgcolor=red\n";
		for (my $j=0;$j<=$#chroms;$j++) {
			open IN, "$outDir/CHG_$chroms[$j].wig";
			while (my $line=<IN>) {print WIGOUT "$line"}
			close IN;
			unlink ("$outDir/CHG_$chroms[$j].wig");
		}
		close WIGOUT;
	}
	if ($context eq "CHH" or $context eq "ALL") {
		open WIGOUT, ">$outDir/CHH.wig";
		print WIGOUT "track type=wiggle_0 name=CHH Methylation  description=MethylExtract_1.9 glyph=wiggle_density bgcolor=red fgcolor=red\n";
		for (my $j=0;$j<=$#chroms;$j++) {
			open IN, "$outDir/CHH_$chroms[$j].wig";
			while (my $line=<IN>) {print WIGOUT "$line"}
			close IN;
			unlink ("$outDir/CHH_$chroms[$j].wig");
		}
		close WIGOUT;
	}
}

########### Printing LOG FILE ################

print "\nPrinting LOG file\n";

my $seconds_elapsed=time - $time_start;
my $time=&userTime($seconds_elapsed);
open LOG, ">$outDir/Ratios$context"."Stats.log";
print LOG "##########   MethylExtract_1.9 LogFile   ##########\n";
print LOG "\nWORKING FEATURES:\n";
print LOG "Input Directory: $inDir\n";
print LOG "$typeFasta-Fasta Sequences: $seq\n";
print LOG "Output Directory: $outDir\n";
print LOG "Elapsed time: $time\tNumber of threads: $p\tNumber of chromosomes divisions: $chromDiv\n";
print LOG "Strand tags used: Watson strand $tagWout & Crick strand $tagCout\n";
print LOG "Qscore used: $qscore\n";
print LOG "Output: context=$context\tWIG=$wigOut\tBED=$bedOut\n";
print LOG "\nQUALITY FEATURES:\n";
print LOG "Duplicated reads deletion: $delDup\tSimilar nucleotides to consider a duplicated read: $simDupPb\n";
print LOG "First number of positions ignored: $FirstIgnor\n";
print LOG "Last number of positions ignored: $LastIgnor\n";
print LOG "Pair-end overlaping correction: $peOverlap\n";
print LOG "Quality filters:\n  minQ=$Q\tminDepthMeth=$minDepthMeth\tminDepthSNV=$minDepthSNV\n  methNonCpGs=$methNonCpGs\tmaxStrandBias=$maxStrandBias\tvarFraction=$minSNVperc\tmaxPval=$maxPval\n";
print LOG "\nQUALITY RESULTS:\n";
print LOG "Deleted duplicated reads: $duplicatedReadsDel\n";
print LOG "Discarded reads by bisulfite check: $discardReads\n";
print LOG "Discarded positions Q<$Q: $discardPos\n";
print LOG "Discarded methylation ratios with depth<$minDepthMeth: $discardDepthPos\n";
print LOG "Unchecked variation in positions with depth<$minDepthSNV: $uncheckedDepthSNV\n";
print LOG "Discarded SNVs with strand bias over $maxStrandBias: $strandBiasDiscards\n";
print LOG "\nSNVs RESULTS:\n";
foreach (@chroms) {
	print LOG " SNVs on $_:\n";
	if ($homoSNV{$_}) {print LOG "  Homozygous SNVs: $homoSNV{$_}\n"}
	else {print LOG "  Homozygous SNVs: 0\n"}
	if ($heteroSNV{$_}) {print LOG "  Heterozygous SNVs: $heteroSNV{$_}\n"}
	else {print LOG "  Heterozygous SNVs: 0\n"}
}
print LOG "\nMethylation RESULTS:\n";
if ($context eq 'CG') {
	&outMethResults("CG",$coverContextCG_ref,\%numContext,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
}
elsif ($context eq 'CHG') {
	&outMethResults("CHG",$coverContextCHG_ref,\%numContext,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
}
elsif ($context eq 'CHH') {
	&outMethResults("CHH",$coverContextCHH_ref,\%numContext,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
}
elsif ($context eq 'ALL') {
	&outMethResults("CG",$coverContextCG_ref,\%numContext,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
	&outMethResults("CHG",$coverContextCHG_ref,\%numContext,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
	&outMethResults("CHH",$coverContextCHH_ref,\%numContext,$countMeth_ref,$countUnmeth_ref,$sumTotal_ref);
}
else {}

close LOG;
unlink("$inDir/" . 'MultiFASTA.fa');

print "\n";
print "# Methylation Ratios & SNVs retrieved, consult your results in $outDir #\n";
print "# See the log file for further information #\n";


#################    SUBPROCESS     #################

#Checking selected options
#0->hash selected contexts 1->Bisulfite checking 2->nonCpG int/fraction 3->qscore type 4->minimum Q
#5->variation fraction 6->output directory 7->BED format 8->WIG format 9->threads 10->duplicate deletion
#11->ignored bases 12->watson FLAG 13->crick FLAG 14->VCF format 15->similar bases
#16->chromosome divisions 17->maximum reads kept in memory 18->sort option 19->minimum coverage 20->pvalue threshold
#21->minimum methylation depth 22->maximum SNV strand bias 23->fasta input type 24->multi fasta IDs 25->peOverlap
#26->ignored last bases
sub GetOptions{
my %opts;
#Checking General Arguments & Help
if (@ARGV) {
	foreach (@ARGV) {
		my @GetOpt=split(/=/,$_);
		if ($GetOpt[0] eq "flagW" or $GetOpt[0] eq "flagC") {
			if ($GetOpt[1] eq "0") {$opts{$GetOpt[0]}="X"}
			else {$opts{$GetOpt[0]}=$GetOpt[1]}
		}
		elsif ($GetOpt[0] eq "tagW" or $GetOpt[0] eq "tagC") { 	# Deprecated
			if ($GetOpt[1] eq "0") {$opts{$GetOpt[0]}="X"}		# Deprecated
			else {$opts{$GetOpt[0]}=$GetOpt[1]}					# Deprecated
		}														# Deprecated
		else {$opts{$GetOpt[0]}=$GetOpt[1]}
	}
	#Checking input options
	foreach my $keyOpts (keys %opts){
		if (!exists($optionsDefault{$keyOpts})) {die "$keyOpts is not an accepted parameter, please check spelling and case sensitive\n";}
		else {}
	}
	#Checking sequences and indexed files paths
	if($opts{inDir}){
		$inDir=$opts{inDir};
		if (-d $inDir) {
			#Checking sam files
			my $checkAlign="N";
			foreach (glob "$inDir/*.*") {
				if (-e $_) {
					$_ =~ m/(\w+)$/;
					if ($1 eq "sam" or $1 eq "bam" or $1 eq "gz") {$checkAlign="Y";}
					else {}
				}
				else {}
			}
			if ($checkAlign eq "N") {die "Don't exist *.sam, *.bam or *.gz files in the directory\n";}
			else {}
		}
		else {die "Cannot find Alignments Directory: $inDir\n";}
	}
	else {die "Use inDir=[alignments directory] to specify alignment path\n";}
	#Getting patterns to extract
	if ($opts{context}){
		$_[0]=$opts{context};
		if (uc($opts{context}) eq "ALL") {}
		elsif (uc($opts{context}) eq "CG") {}
		elsif (uc($opts{context}) eq "CHG") {}
		elsif (uc($opts{context}) eq "CHH") {}
		else {die "Unavailable context option: use context=[context] with CG,CHG,CHH or ALL\n";}
	}
	else {$_[0]="CG"}
	#Delete duplicates
	if ($opts{delDup}) {
		if (uc($opts{delDup}) eq "N") {$_[10]=uc($opts{delDup});}
		elsif (uc($opts{delDup}) eq "Y") {$_[10]=uc($opts{delDup});}
		else {die "Unavailable delDup option: use delDup= Y or N\n";}
	}
	else {$_[10]="N";}
	#Checking tags
	my (%tagsW,%tagsC);
	if ($opts{flagW}) {
		if ($opts{flagW} eq "X") {$tagsW{0}=1}
		else {%tagsW  = map { $_ => 1 } split(/,/, $opts{flagW})}
	}
	elsif ($opts{tagW}) {																	#	Deprecated
		if ($opts{tagW} eq "X") {$tagsW{0}=1}												#	Deprecated
		else {%tagsW  = map { $_ => 1 } split(/,/, $opts{tagW})}							#	Deprecated
		print "\nWarning: tagW is a deprecated parameter, it's been replaced by flagW\n";	#	Deprecated
	}																						#	Deprecated
	else {die "Watson FLAGs are required:\n Common single-end FLAG -> flagW=0\n Common pair-end FLAGs -> flagW=99,147\n"}
	if ($opts{flagC})  {
		if ($opts{flagC} eq "X") {$tagsC{0}=1}
		else {%tagsC = map { $_ => 1 } split(/,/, $opts{flagC})}
	}
	elsif ($opts{tagC}) {																	#	Deprecated
		if ($opts{tagC} eq "X") {$tagsW{0}=1}												#	Deprecated
		else {%tagsC  = map { $_ => 1 } split(/,/, $opts{tagC})}							#	Deprecated
		print "\nWarning: tagC is a deprecated parameter, it's been replaced by flagC\n";	#	Deprecated
	}																						#	Deprecated
	else {die "Crick FLAGs are required:\n Common single-end FLAG -> flagC=16\n Common pair-end FLAGs -> flagC=83,163\n"}
	$_[12]=\%tagsW;
	$_[13]=\%tagsC;
	#Number of nonCpG contexts methylated to discard read
	if (defined($opts{methNonCpGs})) {
		if ($opts{methNonCpGs}=~m/[0-9]+/) {
			$_[2]=$opts{methNonCpGs};
			$_[1]="Y";
		}
		elsif ($opts{methNonCpGs}==0) {
			$_[2]=$opts{methNonCpGs};
			$_[1]="N";
		}
		else {die "The number of nonCpG context methylated to discard the read must be a fraction between 0-1 or an integer higher than 1: use methNonCpGs=[real]\n";}
	}
	else {
		$_[1]="Y";
		$_[2]=0.9;
	}
	#Checking Qscore
	if ($opts{qscore}) {
		$qscore=$opts{qscore};
		if ($qscore eq "phred33-quals" or $qscore eq "phred64-quals" or $qscore eq "solexa-quals" or $qscore eq "NA"){}
		else {
			print "Quality score isn't an accepted format, launching by default phred33-quals\n";
			$_[3]="phred33-quals";
		}
	}
	else {$_[3]="phred33-quals";}
	#Minimun base quality
	if (defined($opts{minQ})) {
		$_[4]=$opts{minQ};
		if ($_[4]=~m/[0-9]+/) {
			$_[4]=~s/,/\./;
			if ($_[4]-int($_[4])>0) {die "Minimun base quality must be an integer\n";}
			else {}
		}
		else {die "Minimun base quality must be an integer\n";}
	}
	else {$_[4]=20;}
	#Number of similar bases in a start coordinate similar reads to be considered a duplicated reads
	if (defined($opts{simDupPb})) {
		$_[15]=$opts{simDupPb};
		if ($_[15]=~m/[0-9]+/) {
			$_[15]=~s/,/\./;
			if ($_[15]-int($_[15])>0) {die "Number of similar bases must be an integer\n";}
			else {}
		}
		else {die "Number of similar bases must be an integer\n";}
	}
	else {$_[15]=32;}
	#chromosome divisions
	if (defined($opts{chromDiv})) {
		$_[16]=$opts{chromDiv};
		if ($_[16]=~m/[0-9]+/) {
			$_[16]=~s/,/\./;
			if ($_[16]-int($_[16])>0) {die "Chromosomes divisions must be an integer\n";}
			else {}
		}
		else {die "Chromosomes divisions must be an integer\n";}
		if ($_[16]==0) {die "Chromosomes divisions cannot be an 0\n";}
		else {}
	}
	else {$_[16]=400;}
	#number of temp lines to keep in memory
	if (defined($opts{memNumReads})) {
		$_[17]=$opts{memNumReads};
		if ($_[17]=~m/[0-9]+/) {
			$_[17]=~s/,/\./;
			if ($_[17]-int($_[17])>0) {die "Number of temporally reads to keep in memory must be an integer\n";}
			else {}
		}
		else {die "Number of temporally reads to keep in memory must be an integer\n";}
		if ($_[17]==0) {die "Number of temporally reads to keep in memory cannot be an 0\n";}
		else {}
	}
	else {$_[17]=200000;}
	#First number of positions ignored
	if (defined($opts{FirstIgnor})) {
		if ($opts{FirstIgnor}==0) {$_[11]=0}
		else {$_[11]=$opts{FirstIgnor}}
		if ($_[11]=~m/[0-9]+/) {
			$_[11]=~s/,/\./;
			if ($_[11]-int($_[11])>0) {die "Number of first based ignored must be an integer\n";}
			else {}
		}
		else {die "Number of first based ignored must be a real integer\n";}
	}
	else {$_[11]=0;}
	#Last number of positions ignored
	if (defined($opts{LastIgnor})) {
		if ($opts{LastIgnor}==0) {$_[26]=0}
		else {$_[26]=$opts{LastIgnor}}
		if ($_[26]=~m/[0-9]+/) {
			$_[26]=~s/,/\./;
			if ($_[26]-int($_[26])>0) {die "Number of last based ignored must be an integer\n";}
			else {}
		}
		else {die "Number of last based ignored must be a real integer\n";}
	}
	else {$_[26]=0;}
	#Minimum depth SNVs
	if (defined($opts{minDepthSNV})) {
		$_[19]=$opts{minDepthSNV};
		if ($_[19]=~m/[0-9]+/) {
			$_[19]=~s/,/\./;
			if ($_[19]-int($_[19])>0) {die "Minimum depth for SNVs calls must be an integer\n";}
			else {}
		}
		else {die "Minimum depth for SNV calls must be an integer\n";}
		if ($_[19]<=0) {die "Minimum depth for SNV calls must be higher than 0\n";}
	}
	else {$_[19]=1;}
	#Minimum depth Methylation
	if (defined($opts{minDepthMeth})) {
		$_[21]=$opts{minDepthMeth};
		if ($_[21]=~m/[0-9]+/) {
			$_[21]=~s/,/\./;
			if ($_[21]-int($_[21])>0) {die "Minimum depth for Methylation calls must be an integer\n";}
			else {}
		}
		else {die "Minimum depth for Methylation calls must be an integer\n";}
		if ($_[21]<=0) {die "Minimum depth for Methylation calls must be higher than 0\n";}
	}
	else {$_[21]=1;}
	#Perc of SNVs to discard positions
	if (defined($opts{varFraction})) {
		$_[5]=$opts{varFraction};
		if ($_[5]=~m/[0-9]+/) {
			$_[5]=~s/,/\./;
			if (int($_[5])>1) {die "Minimum allele frequency must be between 0-1\n";}
			elsif ($_[5]==1) {}
			else {}
		}
		else {die "Minimum allele frequency must be between 0-1\n";}
	}
	else {$_[5]=0.1;}
	#Maximum allowed strand bias
	if (defined($opts{maxStrandBias})) {
		$_[22]=$opts{maxStrandBias};
		if ($_[22]=~m/[0-9]+/) {
			$_[22]=~s/,/\./;
			if (int($_[22])>1) {die "Maximum strand bias must be between 0-1\n";}
			elsif ($_[22]==1) {}
			else {}
		}
		else {die "Maximum strand bias must be between 0-1\n";}
	}
	else {$_[22]=0.7;}
	#Checking output directory
	if ($opts{outDir}) {
		$_[6]=$opts{outDir};
		if (-d $_[6]) {}
		else {mkdir($outDir);}
	}
	else {$_[6]=$inDir;}
	#BED output
	if ($opts{bedOut}) {
		if (uc($opts{bedOut}) eq "N") {
			$_[7]=uc($opts{bedOut});
		}
		elsif (uc($opts{bedOut}) eq "Y") {
			$_[7]=uc($opts{bedOut});
		}
		else {die "Unavailable bedOut option: use bedOut= Y or N\n";}
	}
	else {$_[7]="N";}
	#WIG output
	if ($opts{wigOut}) {
		if (uc($opts{wigOut}) eq "N") {
			$_[8]=uc($opts{wigOut});
		}
		elsif (uc($opts{wigOut}) eq "Y") {
			$_[8]=uc($opts{wigOut});
		}
		else {die "Unavailable bedOut option: use wigOut= Y or N\n";}
	}
	else {$_[8]="N";}
	#VCF output
	if ($opts{vcfOut}) {
		if (uc($opts{vcfOut}) eq "N") {
			$_[14]=uc($opts{vcfOut});
		}
		elsif (uc($opts{vcfOut}) eq "Y") {
			$_[14]=uc($opts{vcfOut});
		}
		else {die "Unavailable vcfOut option: use vcfOut= Y or N\n";}
	}
	else {$_[14]="N";}
	#input sam files splitted by chromosome
	if ($opts{chromSplitted}) {
		if (uc($opts{chromSplitted}) eq "N") {
			$_[18]=uc($opts{chromSplitted});
		}
		elsif (uc($opts{chromSplitted}) eq "Y") {
			$_[18]=uc($opts{chromSplitted});
		}
		else {die "Unavailable chromSplitted option: use chromSplitted= Y or N\n";}
	}
	else {$_[18]="N";}
	#input sam files splitted by chromosome																	#	Deprecated
	if ($opts{chromSorted}) {																				#	Deprecated
		print "\nWarning: chromSorted is a deprecated parameter, it's been replaced by chromSplitted\n";	#	Deprecated
		if (uc($opts{chromSorted}) eq "N") {																#	Deprecated
			$_[18]=uc($opts{chromSorted});																	#	Deprecated
		}																									#	Deprecated
		elsif (uc($opts{chromSorted}) eq "Y") {																#	Deprecated
			$_[18]=uc($opts{chromSorted});																	#	Deprecated
		}																									#	Deprecated
		else {die "Unavailable chromSplitted option: use chromSplitted= Y or N\n";}							#	Deprecated
	}																										#	Deprecated
	else {}																									#	Deprecated
	#Checking threads number
	if ($opts{p}) {
		$_[9]=$opts{p};
		if ($_[9]=~m/[0-9]+/) {
			$_[9]=~s/,/\./;
			if ($_[9]-int($_[9])>0) {die "Number of Pipe threads must be an integer\n";}
			elsif ($_[9]<1) {die "Number of Pipe threads must be 1 or more\n";}
			else {}
		}
		else {die "Number of threads must be an integer\n";}
	}
	else {$_[9]=4;}
	#Pvalue threshold
	if (defined($opts{maxPval})) {
		$_[20]=$opts{maxPval};
		if ($_[20]=~m/[0-9]+/) {
			$_[20]=~s/,/\./;
			if (int($_[20])>1) {die "Pvalue threshold must be between 0 and 1\n";}
			elsif ($_[20]==1) {}
			else {}
		}
		else {die "Pvalue threshold must be between 0 and 1\n";}
	}
	else {$_[20]=0.05;}
	#Checking sequence path
	if($opts{seq}){
		$seq=$opts{seq};
		if (-d $seq) {
			#Checking fasta sequences
			foreach (glob "$seq/*.{fa,fas,fasta,fna}") {
				if (-e $_) {
					open testFA, "$_";
					my $line=<testFA>;
					if (substr($line,0,1) eq ">") {}
					else {die "$_ doesn't seem to be a fasta file\n"}
					close testFA;
				}
				else {die "Doesn't exist *.fa, *.fas, *.fasta or *.fna files on the directory\n";}
			}
			$_[23]="Single";
		}
		elsif (-f $seq) {
			#checking multifasta
			print "\nChecking and charging Multi-fasta IDs...\n";
			
			open (INFILE, "< $seq") or die "Can't open file";
			while (<INFILE>) {
				my $line = $_;
				chomp $line;
				if ($line =~ /\>/) { 
					close OUTFILE;
					my $new_file = "$inDir/" . substr($line,1) . ".fa.tmp";
					open (OUTFILE, ">$new_file") or die "Can't open: $new_file $!";
				}
				print OUTFILE "$line\n";
			}
			close OUTFILE;
			open (MULTI, ">>", "$inDir/" . 'MultiFASTA.fa') or die "Can't open file";
			foreach (glob "$inDir/*.fa.tmp") {
				open (FASTA, "<", $_) or die "Can't open: $_";
				while ( my $line = <FASTA> ) {
					print MULTI $line;
				}
				close FASTA;
				unlink($_);
			}
			close MULTI;
			$seq = "$inDir/" . 'MultiFASTA.fa';
			$_[24]=&GetMultiFastaID($seq);
			if (scalar{$_[24]}>=1) {$_[23]="Multi";}
			else {die "$_ doesn't seem to be a multi-fasta file\n"}
		}
		else {}
	}
	else {die "Use seq=[sequences directory or multifasta single file] to specify genome sequences\n";}
	#PE overlap
	if ($opts{peOverlap}) {
		
		if (uc($opts{peOverlap}) eq "N") {$_[25]=uc($opts{peOverlap});}
		elsif (uc($opts{peOverlap}) eq "Y") {$_[25]=uc($opts{peOverlap});}
		else {die "Unavailable peOverlap option: use peOverlap= Y or N\n";}
	}
	else {$_[25]="N";}
	#Samtools
	if (defined $opts{samtools}) {
		$_[27] = $opts{samtools};
	}
}
else {
	print "\n################   MethylExtract   ###############\n";
	print "##############   Command-line help   #############\n\n";
	print "Launch as:\n  perl MethylExtract.pl seq=<sequences directory or multifasta single file> inDir=<alignments' directory> flagW=<Watson FLAGs (multiple FLAGs comma separated)> flagC=<Crick FLAGs (multiple FLAGs comma separated)> [OPTIONS]\n\n";
	print "Optional Quality parameters:\n";
	print "  qscore=<fastq quality score: phred33-quals, phred64-quals,solexa-quals, solexa1.3-quals or NA> [default: phred33-quals]\n";
	print "  delDup=<delete duplicated reads: Y or N> [default: N]\n  simDupPb=<number of similar nucleotides to detect a duplicated read> [default: 32]\n";
	print "  FirstIgnor=<number of first bases ignored (5' end)> [default: 0]\n  LastIgnor=<number of last bases ignored (3' end)> [default:0]\n";
	print "  peOverlap=<discard second mate overlapping segment on pair-end alignment reads: Y or N> [default: N]\n";
	print "  minDepthMeth=<minimum number of reads requiered to consider a methylation value in a certain position> [default: 1]\n  minDepthSNV=<minimum number of reads requiered to consider a SNV value in a certain position> [default: 1]\n";
	print "  minQ=<minimun PHRED quality per sequenced nucleotide> [default: 20]\n";
	print "  methNonCpGs=<nonCpG contexts methylated to discard read> [default: 0.9] (methNonCpGs=0 deactivates bisulfite read check)\n  varFraction=<Minimum allele frequency> [default: 0.1]\n";
	print "  maxStrandBias=<Maximum strand bias> [default: 0.7] (maxStrandBias=0 deactivates the threshold)\n";
	print "  maxPval=<Variation p-value threshold> [default: 0.05]\n";
	print "Optional working parameters:\n";
	print "  p=<threads number> [default: 4]\n  chromDiv=<number of chromosome divisions to sort reads> [default: 400]\n";
	print "  memNumReads=<number of lines kept on memory for each thread> [default: 200000]\n";
	print "  chromSplitted=<skip alignment chromosome splitting, files must be chromosome splitted and named by chromosome (example: chr1.sam,etc...)> [default: N]\n";
	print "Optional output parameters:\n";
	print "  context=<methylation context to extract: CG, CHG, CHH or ALL> [default: CG]\n";
	print "  outDir=<output directory> [default: inDir]\n  bedOut=<methylation output in BED format> [default: N]\n";
	print "  wigOut=<methylation output in WIG format> [default: N]\n";
	print "  samtools=<path to samtools>\n";
	die "\n";
}
return($inDir,$seq,$_[0],$_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7],$_[8],$_[9],$_[10],$_[11],
	$_[12],$_[13],$_[14],$_[15],$_[16],$_[17],$_[18],$_[19],$_[20],$_[21],$_[22],$_[23],$_[24],
	$_[25],$_[26],$_[27]);
}

#reading sequences
#0->fasta file
sub GetFasta {
	open (SEC, "$_[0]") or die "FASTA files do not seem to have chromosome alignments tags name\n";
  	my $seqst_temp = "";
    my $z = <SEC>;
    my $tes = substr($z,0,1);
    if($tes ne ">"){die "Sequence seems not to be in FASTA format\n";}
    else {}
    $z=~ s/>//;
    my @splitID=split(/\s/,$z);
    my $id=$splitID[0];
    $z =~ s/[\n\t\f\r\s]//g;
    while($z = <SEC>){
      	$z =~ s/[\n\t\f\r_0-9\s]//g;
      	$seqst_temp .= $z;
    }
    return($seqst_temp,$id);
}

#reading sequences
#0->fasta file
sub GetMultiFastaID {
	my @tempIDs=();
	open (SEC, "$_[0]") or die "Multi-FASTA file doesn't exist\n";
	while (my $lineID=<SEC>) {
		chomp($lineID);
		if ($lineID=~s/^>//i) {
			my @splitID=split(/\s/,$lineID);
			my $id=$splitID[0];
			push(@tempIDs,$id);
		}
		else {}
	}
	close SEC;
	return(\@tempIDs);
}

#reading sequences
#0->chromID
sub GetMultiFastaSeq {
	open (SEC, "$seq") or die "Multi-FASTA file doesn't exist\n";
	my $chromCheck=0;
	my $seqst_temp="";
	while (my $lineID=<SEC>) {
		if ($lineID=~/^>$_[0]/i) {
			$chromCheck=1;
			next;
		}
		if ($chromCheck==1) {
			if ($lineID=~/^>/i) {last}
			else {
				chomp($lineID);
				$seqst_temp.=$lineID;
			}
		}
		else {}
	}
	close SEC;
	return($seqst_temp,$_[0]);
}

#Using threads for split datasets
#0->Queue 1->Results 2->Input directory 3->splitChroms function 4->Sequence directory
sub StaticthreadsSplit {
	my $returnData=\my @returnData;
	my %test=();
	if ($typeFasta eq "Single") {
		foreach (glob "$_[4]/*.{fa,fas,fasta,fna}") {
			$_=~/(\S+)\.(fa|fas|fasta|fna)$/i;
			my @splitFile=split(/\//,$1);
			$test{$splitFile[$#splitFile]}=1;
		}
	}
	elsif ($typeFasta eq "Multi") {
		foreach (@{$multiSeqID}) {$test{$_}=1;}
	}
	else {}
	foreach (glob "$_[2]/*.{sam,sam.gz,bam}") {
		$_=~/(\S+)\.(sam|sam.gz|bam)$/i;
		my @splitFile=split(/\//,$1);
		if ($test{$splitFile[$#splitFile]}) {die "Input files seem to be grouped by chromosome, please activate chromSplitted option and launch again\n"}
		$_[0]->enqueue("$_");
	}
	for (1..$p) {
		$_[0]->enqueue(undef);
		push(my @processes, threads->create($_[3]));
	}
	foreach (threads->list){$_->join;}
	$_[1]->enqueue(undef);
	while (my $result = $_[1]->dequeue){
		push(@{$returnData},$result);
	}
	return($returnData);
}

#Using threads for main process
#0->Queue 1->Results 2->Sequence directory 3->Main function
sub Staticthreads {
	my $returnData=\my @returnData;
	if ($typeFasta eq "Single") {
		foreach (glob "$_[2]/*.{fa,fas,fasta,fna}") {$_[0]->enqueue("$_");}
	}
	elsif ($typeFasta eq "Multi") {
		foreach (@{$multiSeqID}) {$_[0]->enqueue("$_");}
	}
	else {}
	for (1..$p) {
		$_[0]->enqueue(undef);
		push(my @processes, threads->create($_[3]));
	}
	foreach (threads->list){$_->join;}
	$_[1]->enqueue(undef);
	while (my $result = $_[1]->dequeue){
    	push(@{$returnData},$result);
	}
	return($returnData);
}

#Calling Qscore from ascii characters
#0->ascii character
sub sanger {return(ord($_[0])-33);}
sub solexa {return((10 * log(1 + 10 ** (ord($_[0]) - 64) / 10)) / log(10));}
sub illumina {return(ord($_[0])-64);}

#Getting read bases
#0->reference sequence 1->SAM line 2->Watson bases hash 3->Crick bases hash
#4->Position bases covered 5->Watson quality hash 6->Crick quality hash 7->reference length
sub getBases {
	my @splitLine=@{$_[1]};
	$splitLine[3]=$splitLine[3]-1;
	my $lenRead=length($splitLine[9]);
	#Validate bisulfite conversion by nonCpG contexts
	my $validate;
	if ($checkBisulfite eq "Y") {
		$validate=&bisulfite(\@splitLine,$lenRead,$_[0],$_[7]);
	}
	else {$validate="Y"}
	#Check Pair-end overlaping ends
	if ($peOverlap eq "Y") {
		if ($splitLine[8]>0 and $splitLine[7]>0) {
			if ($lenRead>$splitLine[8]) {
				$splitLine[9]=substr($splitLine[9],0,$splitLine[8]);
				$splitLine[10]=substr($splitLine[10],0,$splitLine[8]);
				$lenRead=length($splitLine[9]);
			}
			else {}
		}
		else {}
	}
	else {}
	###############################
	if ($validate eq "Y") {
		#Single-End Reads
		if ($splitLine[7]==0 and $splitLine[8]==0) {
			if ($tagW->{$splitLine[1]}) {
				my $ini=$FirstIgnor;
				my $end=$lenRead-($LastIgnor+1);
				($_[2],$_[4],$_[5])=&OrientationgetBases(\@splitLine,$ini,$end,$_[2],$_[3],$_[4],$_[5]);
			}
			elsif ($tagC->{$splitLine[1]}) {
				my $ini=$LastIgnor;
				my $end=$lenRead-($FirstIgnor+1);
				($_[3],$_[4],$_[6])=&OrientationgetBases(\@splitLine,$ini,$end,$_[3],$_[2],$_[4],$_[6]);
			}
			else {}
		}
		#Pair-End Reads
		else {
			if ($tagW->{$splitLine[1]}) {
				#First Mate
				if ($splitLine[8]>0) {
					my $ini=$FirstIgnor;
					my $end=$lenRead-($LastIgnor+1);
					($_[2],$_[4],$_[5])=&OrientationgetBases(\@splitLine,$ini,$end,$_[2],$_[3],$_[4],$_[5]);
				}
				#Second Mate
				else {
					my $ini=$LastIgnor;
					my $end=$lenRead-($FirstIgnor+1);
					($_[2],$_[4],$_[5])=&OrientationgetBases(\@splitLine,$ini,$end,$_[2],$_[3],$_[4],$_[5]);
				}
			}
			elsif ($tagC->{$splitLine[1]}) {
				#First Mate
				if ($splitLine[8]>0) {
					my $ini=$FirstIgnor;
					my $end=$lenRead-($LastIgnor+1);
					($_[3],$_[4],$_[6])=&OrientationgetBases(\@splitLine,$ini,$end,$_[3],$_[2],$_[4],$_[6]);
				}
				#Second Mate
				else {
					my $ini=$LastIgnor;
					my $end=$lenRead-($FirstIgnor+1);
					($_[3],$_[4],$_[6])=&OrientationgetBases(\@splitLine,$ini,$end,$_[3],$_[2],$_[4],$_[6]);
				}
			}
			else {}
		}
	}
	else {
		{lock($discardReads);
		$discardReads++;;}
	}
	return($_[2],$_[3],$_[4],$_[5],$_[6]);
}

#5'->3' reads getBases or 3'->5' reads getBases
#0->\@splitLine 1->$ini 2->$end
#3->Strand bases hash 4->RevStrand bases hash 
#5->Position bases covered 6->Strand Quality hash
sub OrientationgetBases {
	#Reading CIGAR
	my $CIGAR="";
	my @count = ${$_[0]}[5] =~ /(\d+)/g;
	my @chars = ${$_[0]}[5] =~ /(\D+)/g;
	for (my $i=0;$i<=$#count;$i++) {$CIGAR.=$chars[$i]x$count[$i]}
	my $CIGARpos=${$_[0]}[3];
	for (my $i=0;$i<=length($CIGAR)-1;$i++) {
		if (substr($CIGAR,$i,1)=~/[M=X]/) {
			my $base=uc(substr(${$_[0]}[9],0,1));
			substr(${$_[0]}[9],0,1)='';
			my $ascii=substr(${$_[0]}[10],0,1);
			substr(${$_[0]}[10],0,1)='';
			my $PHREDval=&{$subQscore{$qscore}}($ascii);
			if ($PHREDval>=$Q) {($_[3],$_[5],$_[6])=&store($CIGARpos,$_[3],$base,$_[4],$_[5],$_[6],$ascii);}
			else {
				{lock($discardPos);
				$discardPos++;}
			}
			$CIGARpos++;
		}
		elsif (substr($CIGAR,$i,1)=~/[DN]/) {$CIGARpos++}
		elsif (substr($CIGAR,$i,1)=~/[PH]/) {}
		elsif (substr($CIGAR,$i,1)=~/[IS]/) {
			substr(${$_[0]}[9],0,1)='';
			substr(${$_[0]}[10],0,1)='';
		}
		else {
			print "${$_[0]}[5] skipped, unknown CIGAR character\n";
			last;
		}
	}
	return($_[3],$_[5],$_[6]);
}

#Checking bisulfite conversion
#0->\@splitLine 1->read length 2->reference seq 3->reference length
sub bisulfite {
	if ($methNonCpGs>0) {
		my $NoCpGseq=0;
		my $methNoCpGread=0;
		#Single-End
		if (${$_[0]}[7]==0 and ${$_[0]}[8]==0) {
			if ($tagW->{${$_[0]}[1]}) {
				my $ini=$FirstIgnor;
				my $end=$_[1]-$LastIgnor-1;
				($NoCpGseq,$methNoCpGread)=&OrientationBisulfite($_[0],$_[2],$_[1],$NoCpGseq,$methNoCpGread,"C","G",$ini,$end,1,$_[3]);
			}
			elsif ($tagC->{${$_[0]}[1]}) {
				my $ini=$LastIgnor;
				my $end=$_[1]-$FirstIgnor-1;
				($NoCpGseq,$methNoCpGread)=&OrientationBisulfite($_[0],$_[2],$_[1],$NoCpGseq,$methNoCpGread,"G","C",$ini,$end,-1,$_[3]);
			}
			else {return("N")}
		}
		#Pair-End
		else {
			if ($tagW->{${$_[0]}[1]}) {
				if (${$_[0]}[8]>0) {
					my $ini=$FirstIgnor;
					my $end=$_[1]-$LastIgnor-1;
					($NoCpGseq,$methNoCpGread)=&OrientationBisulfite($_[0],$_[2],$_[1],$NoCpGseq,$methNoCpGread,"C","G",$ini,$end,1,$_[3]);
				}
				else {
					my $ini=$LastIgnor;
					my $end=$_[1]-$FirstIgnor-1;
					($NoCpGseq,$methNoCpGread)=&OrientationBisulfite($_[0],$_[2],$_[1],$NoCpGseq,$methNoCpGread,"C","G",$ini,$end,1,$_[3]);
				}
			}
			elsif ($tagC->{${$_[0]}[1]}) {
				if (${$_[0]}[8]>0) {
					my $ini=$FirstIgnor;
					my $end=$_[1]-$LastIgnor-1;
					($NoCpGseq,$methNoCpGread)=&OrientationBisulfite($_[0],$_[2],$_[1],$NoCpGseq,$methNoCpGread,"G","C",$ini,$end,-1,$_[3]);
				}
				else {
					my $ini=$LastIgnor;
					my $end=$_[1]-$FirstIgnor-1;
					($NoCpGseq,$methNoCpGread)=&OrientationBisulfite($_[0],$_[2],$_[1],$NoCpGseq,$methNoCpGread,"G","C",$ini,$end,-1,$_[3]);
				}
			}
			else {return("N")}
		}
		if ($methNonCpGs>0 and $methNonCpGs<=1) {
			my $percMeth=0;
			if ($NoCpGseq>0) {$percMeth=$methNoCpGread/$NoCpGseq;}
			else {}
			if ($percMeth<$methNonCpGs) {return("Y");}
			else {return("N")}
		}
		elsif ($methNonCpGs>1) {
			if ($methNoCpGread>$methNonCpGs) {return("N");}
			else {return("Y")}
		}
		else {}
	}
	else {return("Y")}
}

#5'->'3 bisulfite 3'->5' bisulfite
#0->SAM line 1->reference sequence, 2->read length, 3->NoCpGseq 4->methNoCpGread 
#5->methBase 6->CpG (G base) 7->ini 8->$end 9->W(+1) C(-1) 10-> reference length
sub OrientationBisulfite {
	for (my $i=$_[7];$i<=$_[8];$i++) {
		if ($i+${$_[0]}[3]+$_[9]<$_[10] and $i+${$_[0]}[3]+$_[9]>=0) {
			if (uc(substr($_[1],$i+${$_[0]}[3],1)) eq "$_[5]" and uc(substr($_[1],$i+${$_[0]}[3]+$_[9],1)) ne "$_[6]") {
				$_[3]++;
				if (uc(substr(${$_[0]}[9],$i,1)) eq "$_[5]") {
					$_[4]++;
				}
				else {}
			}
			else {}
		}
		else {}
	}
	return($_[3],$_[4]);
}

#Store Data
#0->positions 1->hash bases strand (to be kept) 2->sequenced base 3->hash bases strand (to be checked)
#4->array positions 5->hash quals strand (to be kept) 6->ascii character
sub store {
	if ($_[2] ne "N") {
		if ($_[1]->{"$_[0]"}) {
			$_[1]->{"$_[0]"}.="$_[2]";
			$_[5]->{"$_[0]"}.="$_[6]";
		}
		else {
			$_[1]->{"$_[0]"}="$_[2]";
			$_[5]->{"$_[0]"}="$_[6]";
			if ($_[3]->{"$_[0]"}) {}
			else {push(@{$_[4]},$_[0])}
		}
		return($_[1],$_[4],$_[5]);
	}
	else {return($_[1],$_[4],$_[5]);}
}

#Changing time format
sub userTime {
	my $days=int($_[0]/86400);
	$_[0]-=$days*86400;
	my $hours=int($_[0]/3600);
	$_[0]-=$hours*3600;
	my $mins=int($_[0]/60);
	my $segs=$_[0]-($mins*60);
	my $time="$days"."d : $hours"."h : $mins"."m :$segs"."s";
	return($time);
}

#print VCF file
#0-> chromosome 1-> position 2-> reference base 3-> sample bases
#4-> hash basesW 5-> hash quals 6->hash variation 7-> pvalue 8-> hash basesC
sub printVCF {
	my $altBases="";
	my $freqBases="";
	my $qBases="";
	my $depthW=0;
	my $depthC=0;
	my $depthQ=0;
	my $qual=0;
	my $GTbases="";
	my $GTcode="";
	for (my $i=0;$i<=length($_[3])-1;$i++) {
		$depthW+=$_[4]->{substr($_[3],$i,1)};
		$depthC+=$_[8]->{substr($_[3],$i,1)};
		$depthQ+=length($_[5]->{substr($_[3],$i,1)});
		for (my $j=0;$j<=length($_[5]->{substr($_[3],$i,1)})-1;$j++) {$qual+=&{$subQscore{$qscore}}(substr($_[5]->{substr($_[3],$i,1)},$j,1))}
	}
	my $meanQual=sprintf("%.2f",$qual/$depthQ);
	my %alts=();
	my %altQ=();
	for (my $i=0;$i<=length($_[3])-1;$i++) {
		my $alt=substr($_[3],$i,1);
		my $freq=sprintf("%.2f", ($_[4]->{$alt}+$_[8]{$alt})/($depthW+$depthC));
		my $altqual=0;
		for (my $j=0;$j<=length($_[5]->{$alt})-1;$j++) {$altqual+=&{$subQscore{$qscore}}(substr($_[5]->{$alt},$j,1))}
		$alts{$alt}=$freq;
		$altQ{$alt}=sprintf("%.2f", $altqual/length($_[5]->{$alt}));
	}
	my $count=0;
	my $undefGT=0;
	my $countNonRefW=0;
	my $countNonRefC=0;
	for my $keyBases (sort {$alts{$b} <=> $alts{$a}} keys %alts) {
		if ($keyBases ne $_[2]) {
			$countNonRefW+=$_[4]->{$keyBases};
			$countNonRefC+=$_[8]->{$keyBases};
			$altBases.="$keyBases,";
			$freqBases.="$alts{$keyBases},";
			$qBases.="$altQ{$keyBases},";
			$count++;
		}
		else {$qBases.="$altQ{$keyBases},"."$qBases"}
		if (length($GTbases)<2) {
			$GTbases.=$keyBases;
			if ($keyBases ne $_[2]) {$GTcode.=$count}
		}
		else {if ($alts{substr($GTbases,1,1)}==$alts{$keyBases}) {$undefGT=1}}
	}
	my $countRefW=$depthW-$countNonRefW;
	my $countRefC=$depthC-$countNonRefC;
	$altBases=~s/,$//;
	$freqBases=~s/,$//;
	$qBases=~s/,$//;
	
	my $GT="";
	if (length($GTcode)==1) {
		if (length($GTbases)==1) {$GT="$GTcode"."/"."$GTcode"}
		else {$GT="0"."/"."$GTcode"}
	}
	else {
		my $a=substr($GTcode,0,1);
		my $b=substr($GTcode,1,1);
		$GT="$a"."/"."$b";
	}
	if ($minSNVperc<=0.5) {open VCF, ">>$outDir/SNVs_$_[0].output"}
	my $pos1=$_[1]+1;
	my $QUAL=1000;
	if ($_[7]>0) {$QUAL=-10*log($_[7]);}
	else {}
	$QUAL=sprintf("%.2f",$QUAL);
	my $depth=$depthW+$depthC;
	my $SB=2;
	if ($countRefW==0 and $countNonRefW==0) {}
	elsif ($countRefC==0 and $countNonRefC==0) {}
	else {
		$SB=&TestFisher($countRefW,$countRefC,$countNonRefW,$countNonRefC);
		$SB=1-sprintf("%.3f",$SB);
	}
	if ($undefGT==1) {}
	elsif ($maxStrandBias!=0 and $SB>$maxStrandBias) {
		{lock($strandBiasDiscards);
		$strandBiasDiscards++}
	}
	else {
		if ($minSNVperc<=0.5) {print VCF "$_[0]\t$pos1\t.\t$_[2]\t$altBases\t$QUAL\tPASS\tDP=$depth;AF=$freqBases;MQ=$meanQual;SB=$SB\tGT:BQ:DP4\t$GT:$qBases:$countRefW,$countRefC,$countNonRefW,$countNonRefC\n"}
		$_[6]->{$_[1]}=$GTbases;
		if (length($GTbases)==1) {
			if ($homoSNV{$_[0]}) {
				{lock(%homoSNV);
				$homoSNV{$_[0]}++}
			}
			else {
				{lock(%homoSNV);
				$homoSNV{$_[0]}=1}
			}
		}
		else {
			if ($heteroSNV{$_[0]}) {
				{lock(%heteroSNV);
				$heteroSNV{$_[0]}++}
			}
			else {
				{lock(%heteroSNV);
				$heteroSNV{$_[0]}=1}
			}
		}
	}
	if ($minSNVperc<=0.5) {close VCF}
	return($_[6],$GTbases);
}

#Split Datasets in chromosomes
sub splitChroms {
	while (my $element=$QueueSplit->dequeue()) {
		my $countLines=0;
		my %keptLines=();
		print "Reading $element\n";
		if ($element=~/\.sam$/) {
			open SAM, "$element";
			while (my $line=<SAM>) {
				if ($line!~/^@/i) {
					my @splitLine=split(/\t/,$line);
					$keptLines{$splitLine[2]}{$countLines}=$line;
					$countLines++;
					if ($countLines>=$memNumReads*6) {
						for my $keyChrom (keys %keptLines) {
							open OUT, ">>$inDir/$keyChrom.sam";
							flock(OUT,2);
							for my $keyUni (keys (%{$keptLines{$keyChrom}})) {
								print OUT "$keptLines{$keyChrom}{$keyUni}";
							}
							close OUT;
						}
						%keptLines=();
						$countLines=0;
					}
					else {}
				}
				else {}
			}
			if ($countLines>0) {
				for my $keyChrom (keys %keptLines) {
					open OUT, ">>$inDir/$keyChrom.sam";
					flock(OUT,2);
					for my $keyUni (keys (%{$keptLines{$keyChrom}})) {
						print OUT "$keptLines{$keyChrom}{$keyUni}";
					}
					close OUT;
				}
				%keptLines=();
				$countLines=0;
			}
			else {}
			close SAM;
		}
		elsif ($element=~/\.sam\.gz$/) {
			my $z = new IO::Uncompress::AnyUncompress $element;
			 until (eof($z)) {
			 	my $line = <$z>;
			 	if ($line!~/^@/i) {
			 		my @splitLine=split(/\t/,$line);
			 		$keptLines{$splitLine[2]}{$countLines}=$line;
			 		$countLines++;
			 		if ($countLines>=$memNumReads*6) {
						for my $keyChrom (keys %keptLines) {
							open OUT, ">>$inDir/$keyChrom.sam";
							flock(OUT,2);
							for my $keyUni (keys (%{$keptLines{$keyChrom}})) {
								print OUT "$keptLines{$keyChrom}{$keyUni}";
							}
							close OUT;
						}
						%keptLines=();
						$countLines=0;
					}
					else {}
			 	}
				else {}
			}
			if ($countLines>0) {
				for my $keyChrom (keys %keptLines) {
					open OUT, ">>$inDir/$keyChrom.sam";
					flock(OUT,2);
					for my $keyUni (keys (%{$keptLines{$keyChrom}})) {
						print OUT "$keptLines{$keyChrom}{$keyUni}";
					}
					close OUT;
				}
				%keptLines=();
				$countLines=0;
			}
			else {}
		}
		elsif ($element=~/\.bam$/) {
			#open BAM, "$element";
			open BAM,"$samtools view $element |";
			while (my $line=<BAM>) {
				if ($line!~/^@/i) {
					my @splitLine=split(/\t/,$line);
					$keptLines{$splitLine[2]}{$countLines}=$line;
					$countLines++;
					if ($countLines>=$memNumReads*6) {
						for my $keyChrom (keys %keptLines) {
							open OUT, ">>$inDir/$keyChrom.sam";
							flock(OUT,2);
							for my $keyUni (keys (%{$keptLines{$keyChrom}})) {
								print OUT "$keptLines{$keyChrom}{$keyUni}";
							}
							close OUT;
						}
						%keptLines=();
						$countLines=0;
					}
					else {}
				}
				else {}
			}
			if ($countLines>0) {
				for my $keyChrom (keys %keptLines) {
					open OUT, ">>$inDir/$keyChrom.sam";
					flock(OUT,2);
					for my $keyUni (keys (%{$keptLines{$keyChrom}})) {
						print OUT "$keptLines{$keyChrom}{$keyUni}";
					}
					close OUT;
				}
				%keptLines=();
				$countLines=0;
			}
			else {}
			close BAM;
		}
		else {}
		$ResultsSplit->enqueue();
	}
}

#Main Process
sub Main {
	while (my $element=$Queue->dequeue()) {
		print "Reading $element reference\n";
		my ($seqst,$chrom);
		if ($typeFasta eq "Single") {
			($seqst,$chrom)=&GetFasta("$element");
		}
		elsif ($typeFasta eq "Multi") {
			($seqst,$chrom)=&GetMultiFastaSeq("$element");
		}
		else {}
		my (@UpperPivots,@LowerPivots);
		my $RefLen=length($seqst);
		my $chopSeq=$RefLen/$chromDiv;
		for (my $i=0;$i<=$chromDiv-1;$i++) {
			push(@UpperPivots,$chopSeq*($i+1));
			push(@LowerPivots,$chopSeq*$i);
		}
		print "Extracting alignments from $chrom\n";
		my $maxLenRead=0;
		my ($tempLines,$countLines);
		my $continue="N";
		
		#extracting chromsome alignments
		if (-e "$inDir/$chrom.sam") {
			$continue="Y";
			open CHROM, "$inDir/$chrom.sam";
			$countLines=0;
			while (my $line=<CHROM>) {
				if ($line!~/^@/i) {
					my @splitLine=split(/\t/,$line);
					my $len=length($splitLine[9]);
					if ($len>$maxLenRead) {$maxLenRead=$len};
					if ($delDup eq "N") {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
					else {
						if ($tagW->{$splitLine[1]}) {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
						else {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],$len,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
					}
				}
				else {}
			}
			if ($countLines>0) {
				for (my $i=0;$i<=$#$tempLines;$i++) {
					open TEMP, ">>$inDir/$chrom.$i.sam";
					for (my $j=0;$j<=$#{$$tempLines[$i]};$j++) {
						print TEMP "$$tempLines[$i][$j]";
					}
					close TEMP;
				}
				undef $tempLines;
				$countLines=0;
			}
			else {}
			close CHROM;
			if ($chromSort eq "Y") {}
			else {unlink "$inDir/$chrom.sam"}
		}
		elsif (-e "$inDir/$chrom.sam.gz") {
			my $z = new IO::Uncompress::AnyUncompress "$inDir/$chrom.sam.gz";
			$continue="Y";
				$countLines=0;
				until (eof($z)) {
					my $line = <$z>;
					if ($line!~/^@/i) {
						my @splitLine=split(/\t/,$line);
						my $len=length($splitLine[9]);
						if ($len>$maxLenRead) {$maxLenRead=$len};
						if ($delDup eq "N") {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
						else {
							if ($tagW->{$splitLine[1]}) {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
							else {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],$len,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
						}
					}
					else {}
				}
				if ($countLines>0) {
					for (my $i=0;$i<=$#$tempLines;$i++) {
						open TEMP, ">>$inDir/$chrom.$i.sam";
						for (my $j=0;$j<=$#{$$tempLines[$i]};$j++) {
							print TEMP "$$tempLines[$i][$j]";
						}
						close TEMP;
					}
					undef $tempLines;
					$countLines=0;
				}
				else {}
		}
		elsif (-e "$inDir/$chrom.bam") {
			$continue="Y";
			open BAM,"$samtools view $inDir/$chrom.bam |";
			$countLines=0;
			while (my $line=<BAM>) {
				if ($line!~/^@/i) {
					my @splitLine=split(/\t/,$line);
					my $len=length($splitLine[9]);
					if ($len>$maxLenRead) {$maxLenRead=$len};
					if ($delDup eq "N") {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
					else {
						if ($tagW->{$splitLine[1]}) {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
						else {($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],$len,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);}
					}
				}
				else {}
			}
			if ($countLines>0) {
				for (my $i=0;$i<=$#$tempLines;$i++) {
					open TEMP, ">>$inDir/$chrom.$i.sam";
					for (my $j=0;$j<=$#{$$tempLines[$i]};$j++) {
						print TEMP "$$tempLines[$i][$j]";
					}
					close TEMP;
				}
				undef $tempLines;
				$countLines=0;
			}
			else {}
			close CHROM;
			if ($chromSort eq "Y") {}
			else {unlink "$inDir/$chrom.bam"}
		}
		else {print "It does not seem to exist alignments for $chrom\n"}
		
		if ($continue eq "Y") {
			#Sorting temp files & deleting duplicated reads
			if ($delDup eq "N") {
				print "Sorting $chrom temp files\n";
				open SORT, ">$inDir/$chrom.sort";
				&sortFiles($chrom,\@UpperPivots);
				close SORT;	
			}
			else {
				print "Deleting duplicated reads from $chrom\n";
				open DELDUP, ">$inDir/$chrom.deldup";
				&deleteDup($chrom,\@UpperPivots);
				close DELDUP;
				print "Sorting $chrom temp files\n";
				open NODUP, "$inDir/$chrom.deldup";
				$countLines=0;
				while (my $line=<NODUP>) {
					my @splitLine=split(/\t/,$line);
					($tempLines,$countLines)=&tempFiles($chrom,$splitLine[3],0,\@UpperPivots,\@LowerPivots,$line,$tempLines,$countLines);
				}
				if ($countLines>0) {
					for (my $i=0;$i<=$#$tempLines;$i++) {
						open TEMP, ">>$inDir/$chrom.$i.sam";
						for (my $j=0;$j<=$#{$$tempLines[$i]};$j++) {
							print TEMP "$$tempLines[$i][$j]";
						}
						close TEMP;
					}
					undef $tempLines;
					$countLines=0;
				}
				else {}
				close NODUP;
				unlink("$inDir/$chrom.deldup");
				open SORT, ">$inDir/$chrom.sort";
				&sortFiles($chrom,\@UpperPivots);
				close SORT;	
			}
			
			{lock(@chroms);
			push(@chroms,$chrom);}
			
			#Counting reference contexts
			print "Retrieving methylation contexts in $chrom reference\n";
			if ($context eq "CG" or $context eq "ALL") {
				for (my $i=0;$i<=length($seqst)-2;$i++) {
					if ($CGseq{uc(substr($seqst,$i,2))}) {
						if ($numContext{"$chrom"."_CG"}) {
							{lock(%numContext);
							$numContext{"$chrom"."_CG"}++}
						}
						else {
							{lock(%numContext);
							$numContext{"$chrom"."_CG"}=1}
						}
					} 
				}
			}
			if ($context eq "CHG" or $context eq "ALL") {
				for (my $i=0;$i<=length($seqst)-3;$i++) {
					if ($CHGseq{uc(substr($seqst,$i,3))}) {
						if ($numContext{"$chrom"."_CHG"}) {
							{lock(%numContext);
							$numContext{"$chrom"."_CHG"}++}
						}
						else {
							{lock(%numContext);
							$numContext{"$chrom"."_CHG"}=1}
						}
					} 
				}
			}
			if ($context eq "CHH" or $context eq "ALL") {
				for (my $i=0;$i<=length($seqst)-3;$i++) {
					if ($CHHseq{uc(substr($seqst,$i,3))}) {
						if ($numContext{"$chrom"."_CHH"}) {
							{lock(%numContext);
							$numContext{"$chrom"."_CHH"}++}
						}
						else {
							{lock(%numContext);
							$numContext{"$chrom"."_CHH"}=1}
						}
					} 
				}
			}

			my ($basesW_ref,$basesC_ref,$qualW_ref,$qualC_ref,$basesPos_ref,$cutReads_ref)=();
			print  "Retrieving methylation and SNVs from $chrom\n";
			open READS, "$inDir/$chrom.sort";
				if ($wigOut eq "Y") {
					if ($context eq "CG" or $context eq "ALL") {
						open WIGOUT, ">>$outDir/CG_$chrom.wig";
						print WIGOUT "variableStep\tchrom=$chrom\tspan=2\n";
						close WIGOUT;
					}
					if ($context eq "CHG" or $context eq "ALL") {
						open WIGOUT, ">>$outDir/CHG_$chrom.wig";
						print WIGOUT "variableStep\tchrom=$chrom\tspan=3\n";
						close WIGOUT;
					}
					if ($context eq "CHH" or $context eq "ALL") {
						open WIGOUT, ">>$outDir/CHH_$chrom.wig";
						print WIGOUT "variableStep\tchrom=$chrom\tspan=3\n";
						close WIGOUT;
					}
				}
				my $countLines=0;
				my $endPos=0;
				
				while (my $line=<READS>) {
					my ($nextLine,$nextIni,$actualIni);
					if ($#$cutReads_ref>=0) {
						foreach (@{$cutReads_ref}) {
							my @splitChunkReads=split(/\t/,$_);
							($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref)=&getBases($seqst,\@splitChunkReads,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$RefLen);
						}
						$countLines++;
						$cutReads_ref=();
					}
					else {}
					my @splitLine=split(/\t/,$line);
					($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref)=&getBases($seqst,\@splitLine,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$RefLen);
					my $ReadEnd=$splitLine[3]+(length($splitLine[9])-1);
					if ($ReadEnd>$endPos) {$endPos=$ReadEnd}
					$countLines++;
					if (eof(READS)) {
						&extractData($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$seqst,$chrom);
						$basesW_ref=();
						$basesC_ref=();
						$qualW_ref=();
						$qualC_ref=();
						$basesPos_ref=();
						$countLines=0;
					}
					else {
						if ($countLines==$memNumReads) {
							$actualIni=tell(READS);
							$nextLine=<READS>;
							$nextIni=tell(READS);
							chomp($nextLine);
							my @splitNextLine=split(/\t/,$nextLine);
							if ($endPos>=$splitNextLine[3]) {
								($cutReads_ref,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref)=&chunkReads(\@splitNextLine,$endPos,$seqst,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$cutReads_ref,$nextLine,$RefLen);
								while (1) {
									$actualIni=tell(READS);
									$nextLine=<READS>;
									$nextIni=tell(READS);
									chomp($nextLine);
									@splitNextLine=split(/\t/,$nextLine);
									if (eof(READS)) {
										if ($#$cutReads_ref>=0) {
											foreach (@{$cutReads_ref}) {
												my @splitChunkReads=split(/\t/,$_);
												($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref)=&getBases($seqst,\@splitChunkReads,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$RefLen);
											}
											$countLines++;
											$cutReads_ref=();
										}
										else {}
										($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref)=&getBases($seqst,\@splitNextLine,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$RefLen);
										&extractData($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$seqst,$chrom);
										$basesW_ref=();
										$basesC_ref=();
										$qualW_ref=();
										$qualC_ref=();
										$basesPos_ref=();
										$countLines=0;
										last;
									}
									elsif ($endPos>=$splitNextLine[3]) {
										($cutReads_ref,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref)=&chunkReads(\@splitNextLine,$endPos,$seqst,$basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$cutReads_ref,$nextLine,$RefLen);
									}
									else {
										&extractData($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$seqst,$chrom);
										$basesW_ref=();
										$basesC_ref=();
										$qualW_ref=();
										$qualC_ref=();
										$basesPos_ref=();
										$countLines=0;
										my $backPos=$nextIni-$actualIni;
										seek(READS,-$backPos,1);
										last;
									}
								}
							}
							else {
								&extractData($basesW_ref,$basesC_ref,$basesPos_ref,$qualW_ref,$qualC_ref,$seqst,$chrom);
								$basesW_ref=();
								$basesC_ref=();
								$qualW_ref=();
								$qualC_ref=();
								$basesPos_ref=();
								$countLines=0;
								my $backPos=$nextIni-$actualIni;
								seek(READS,-$backPos,1);
							}
						}
						else {}
					}
				}
			close READS;
		}
		else {}
		unlink "$inDir/$chrom.sort";
		$Results->enqueue();
	}
}

#Writting temp files
#0->chromosome 1->read start position 2->Crick reads correction length 3->Upper pivots array
#4->Lower pivots array 5->Sam line 6->temporally lines 7->temporally lines count
sub tempFiles {
	for (my $i=0;$i<=$#{$_[3]};$i++) {
		if ($_[1]+$_[2]>=$_[4]->[$i] and $_[1]+$_[2]<$_[3]->[$i]) {
			push(@{$_[6][$i]}, $_[5]);
			$_[7]++;
			last;
		}
		else {}
	}
	if ($_[7]>=$memNumReads) {
		for (my $i=0;$i<=$#{$_[6]};$i++) {
			open TEMP, ">>$inDir/$_[0].$i.sam";
			for (my $j=0;$j<=$#{$_[6][$i]};$j++) {
				print TEMP "$_[6][$i][$j]";
			}
			close TEMP;
		}
		$_[6]=();
		$_[7]=0;
	}
	else {}
	return ($_[6],$_[7]);
}

#Sort temp files
#0->chromosome 1->Upper pivots array
sub sortFiles {
	for (my $i=0;$i<=$#{$_[1]};$i++) {
		if (-e "$inDir/$_[0].$i.sam") {
			my %sortFile;
			my $j=0;
			open TEMP, "$inDir/$_[0].$i.sam";
				while (my $line=<TEMP>) {
					my @splitLine=split(/\t/,$line);
					$sortFile{$splitLine[3]}{$j}=$line;
					$j++;
				}
			close TEMP;
			for my $keyStart (sort {$a<=>$b} keys %sortFile) {
				for my $keyUni (keys %{$sortFile{$keyStart}}) {
					print SORT "$sortFile{$keyStart}{$keyUni}";
				}
			}
			%sortFile=();
			unlink ("$inDir/$_[0].$i.sam");
		}
		else {}		
	}	
}

#Deleting duplicated reads
#0->chromosome 1->Upper pivots array
sub deleteDup {
	for (my $i=0;$i<=$#{$_[1]};$i++) {
		if (-e "$inDir/$_[0].$i.sam") {
			my %delDupFile;
			open TEMP, "$inDir/$_[0].$i.sam";
				my $countALLreads=0;
				while (my $line=<TEMP>) {
					my @splitLine=split(/\t/,$line);
					my $len=length($splitLine[9]);
					if ($tagC->{$splitLine[1]}) {$delDupFile{$splitLine[1]}{$splitLine[3]+$len}{$countALLreads}=$line;}
					else {$delDupFile{$splitLine[1]}{$splitLine[3]}{$countALLreads}=$line;}
					$countALLreads++;
				}
			close TEMP;
			for my $keyStrand (sort {$a<=>$b} keys %delDupFile) {
				for my $keyStart (sort {$a<=>$b} keys %{$delDupFile{$keyStrand}}) {
					my @dupReads=();
					for my $keyUni (keys %{$delDupFile{$keyStrand}{$keyStart}}) {push(@dupReads,$delDupFile{$keyStrand}{$keyStart}{$keyUni});}
					if ($#dupReads>0) {
						my %used=();
						for (my $j=0;$j<=$#dupReads;$j++) {
							my @reallyDup=();
							if (!$used{$j}) {
								push(@reallyDup,$dupReads[$j]);
								my @splitArrayA=split(/\t/,$dupReads[$j]);
								$used{$j}=1;
								for (my $k=$j+1;$k<=$#dupReads;$k++) {
									if (!$used{$k}) {
										my @splitArrayB=split(/\t/,$dupReads[$k]);
										if (length($splitArrayA[9])>=length($splitArrayB[9])) {
											my ($dupReads,$reallyDup,$used)=&checkDup($keyStrand,$splitArrayB[9],$splitArrayB[10],$splitArrayA[9],$splitArrayA[10],$k,\@dupReads,\@reallyDup,\%used);
											@dupReads=@{$dupReads};
											@reallyDup=@{$reallyDup};
											%used=%{$used};
										}
										else {
											my ($dupReads,$reallyDup,$used)=&checkDup($keyStrand,$splitArrayA[9],$splitArrayA[10],$splitArrayB[9],$splitArrayB[10],$k,\@dupReads,\@reallyDup,\%used);
											@dupReads=@{$dupReads};
											@reallyDup=@{$reallyDup};
											%used=%{$used};
										}
									}
									else {}
								}
								if ($#reallyDup>0) {
									my %meanQ=();
									for (my $l=0;$l<=$#reallyDup;$l++) {
										my @splitDup=split(/\t/,$reallyDup[$l]);
										my $sumQ=0;
										my $len=length($splitDup[10]);
										for (my $m=0;$m<=length($splitDup[10])-1;$m++) {
											my $ascii=substr($splitDup[10],$m,1);
											my $PHREDval=&{$subQscore{$qscore}}($ascii);
											if ($PHREDval>=$Q) {$sumQ+=1;}
										}
										$meanQ{$sumQ}{$len}=$l;
									}
									for my $q ( sort {$b<=>$a} keys %meanQ) {
										for my $len ( sort {$b<=>$a} keys %{$meanQ{$q}}) {
											print DELDUP "$reallyDup[$meanQ{$q}{$len}]";
											last;
										}
										last;
									}
									%meanQ=();
								}
								else {print DELDUP "$reallyDup[0]"}
								{lock($duplicatedReadsDel);
								$duplicatedReadsDel+=$#reallyDup;}
								@reallyDup=();
							}
							else {}
						}
						%used=();
					}
					else {print DELDUP "$dupReads[0]"}
				}
			}
			%delDupFile=();
			unlink ("$inDir/$_[0].$i.sam");
		}
		else {}		
	}
}

sub checkDup {
	if ($tagW->{$_[0]}) {
		my $seqCheck="";
		my $QCheck="";
		if (length($_[1])>=$simDupPb) {
			$seqCheck=uc(substr($_[1],0,$simDupPb));
			$QCheck=uc(substr($_[2],0,$simDupPb));
		}
		else {
			$seqCheck=uc($_[1]);
			$QCheck=uc($_[2]);
		}
		my $countDiffPB=0;
		for (my $h=0;$h<=length($seqCheck)-1;$h++) {
			my $asciiB=substr($QCheck,$h,1);
			my $PHREDvalB=&{$subQscore{$qscore}}($asciiB);
			my $asciiA=substr($_[4],$h,1);
			my $PHREDvalA=&{$subQscore{$qscore}}($asciiA);
			if ($PHREDvalA>=$Q and $PHREDvalB>=$Q) {
				if (uc(substr($seqCheck,$h,1)) eq uc(substr($_[3],$h,1))) {}
				else {$countDiffPB++}
			}
			else {}
		}
		if ($countDiffPB==0) {
			push(@{$_[7]},${$_[6]}[$_[5]]);
			$_[8]->{"$_[5]"}=1;
		}
		else {}
	}
	else {
		my $seqCheck="";
		my $QCheck="";
		if (length($_[1])>=$simDupPb) {
			$seqCheck=uc(substr($_[1],length($_[1])-$simDupPb,$simDupPb));
			$QCheck=uc(substr($_[2],length($_[2])-$simDupPb,$simDupPb));
		}
		else {
			$seqCheck=uc($_[1]);
			$QCheck=uc($_[2]);
		}
		my $countDiffPB=0;
		for (my $h=1;$h<=length($seqCheck);$h++) {
			my $asciiB=substr($QCheck,length($QCheck)-$h,1);
			my $PHREDvalB=&{$subQscore{$qscore}}($asciiB);
			my $asciiA=substr($_[4],length($_[4])-$h,1);
			my $PHREDvalA=&{$subQscore{$qscore}}($asciiA);
			if ($PHREDvalA>=$Q and $PHREDvalB>=$Q) {
				if (uc(substr($seqCheck,length($seqCheck)-$h,1)) eq uc(substr($_[3],length($_[3])-$h,1))) {}
				else {$countDiffPB++}
			}
			else {}
		}
		if ($countDiffPB==0) {
			push(@{$_[7]},${$_[6]}[$_[5]]);
			$_[8]->{"$_[5]"}=1;
		}
		else {}
	}
	return($_[6],$_[7],$_[8]);
}

sub chunkReads {
	my ($chunkRead,$chunkqRead,$chunkIniRead,$chunkIniqRead,$actualRead,$nextRead);
	my $endChunkRead=$_[0]->[3]+length($_[0]->[9]);
	if ($endChunkRead-1<=$_[1]) {
		($_[3],$_[4],$_[5],$_[6],$_[7])=&getBases($_[2],$_[0],$_[3],$_[4],$_[5],$_[6],$_[7],$_[10]);
	}
	else {
		my $lenActual=($_[1]+1)-$_[0]->[3];
		$chunkIniRead=substr($_[0]->[9],0,$lenActual);
		$chunkIniqRead=substr($_[0]->[10],0,$lenActual);
		my $lenNext=length($_[0]->[9])-$lenActual;
		$chunkRead=substr($_[0]->[9],$_[1]-$_[0]->[3]+1,$lenNext);
		$chunkqRead=substr($_[0]->[10],$_[1]-$_[0]->[3]+1,$lenNext);
		$_[0]->[9]=$chunkIniRead;
		$_[0]->[10]=$chunkIniqRead;
		foreach (@{$_[0]}) {$actualRead.="$_\t"}
		#print TEST "chunkRead:firstPart: $actualRead\n";
		($_[3],$_[4],$_[5],$_[6],$_[7])=&getBases($_[2],$_[0],$_[3],$_[4],$_[5],$_[6],$_[7],$_[10]);
		$_[0]->[3]=$_[1]+1;
		$_[0]->[9]=$chunkRead;
		$_[0]->[10]=$chunkqRead;
		foreach (@{$_[0]}) {$nextRead.="$_\t"}
		chomp $nextRead;
		#print TEST "chunkRead:secondPart: $nextRead\n";
		push(@{$_[8]},$nextRead);
	}
	return($_[8],$_[3],$_[4],$_[5],$_[6],$_[7]);
}

#check SNVs and methylation contexts
#0->watson bases hash 1->crick bases hash 2->positions array 
#3->watson ascii hash 4->crick ascii hash 5->reference sequence
#6->chromosome
sub extractData {
	my $var_ref=();
	my $meth_ref=();
	if (!$_[2]) {
		unlink("$outDir/$_[6].sort");
		print "Chunk without useful sequenced nucleotides on $_[6]\n";
		print "- Please check selected strand FLAGs, Qscore type and minimum Qscore threshold -\n";
		die "\n";
	}
	else {
		foreach (sort {$a<=>$b} @{$_[2]}) {
			my $pos=$_;
			if ($pos>length($_[5])-1) {next}
			my $refBase=uc(substr($_[5],$pos,1));
			my %countBasesW=("A"=>0,"T"=>0,"C"=>0,"G"=>0);
			my %countBasesC=("A"=>0,"T"=>0,"C"=>0,"G"=>0);
			my %countQual=("A"=>"","T"=>"","C"=>"","G"=>"");
			my $countBasesW_ref=\%countBasesW;
			my $countBasesC_ref=\%countBasesC;
			my $countQual_ref=\%countQual;
			my $total=0;
			if ($_[0]->{$pos} and $_[1]->{$pos} and $minSNVperc<1) {
				if (length($_[0]->{$pos})>=$minDepthSNV and length($_[1]->{$pos})>=$minDepthSNV) {
					my $totalBis=0;
					($countBasesW_ref,$countQual_ref,$total)=&countBases($countBasesW_ref,$countQual_ref,$total,$_[0]->{$pos},$_[3]->{$pos},"W");
					($countBasesC_ref,$countQual_ref,$total)=&countBases($countBasesC_ref,$countQual_ref,$total,$_[1]->{$pos},$_[4]->{$pos},"C");
					my $bases="";
					my $countDiffBases=0;
					my $countRefBases=0;
					if ($total>0) {
						($countBasesW_ref,$countBasesC_ref,$totalBis)=&countBasesCorrect($countBasesW_ref,$totalBis,$total,$_[0]->{$pos},$_[3]->{$pos},"W",$countBasesC_ref);
						($countBasesW_ref,$countBasesC_ref,$totalBis)=&countBasesCorrect($countBasesW_ref,$totalBis,$total,$_[1]->{$pos},$_[4]->{$pos},"C",$countBasesC_ref);
						$total+=$totalBis;
						for my $keyBases (keys %{$countBasesW_ref}) {
							if ((($countBasesW_ref->{$keyBases}+$countBasesC_ref->{$keyBases})/$total)>=$minSNVperc) {
								$bases.=$keyBases;
								if ($keyBases ne $refBase) {
									$countDiffBases+=$countBasesW_ref->{$keyBases};
									$countDiffBases+=$countBasesC_ref->{$keyBases};
								}
								else {
									$countRefBases+=$countBasesW_ref->{$keyBases};
									$countRefBases+=$countBasesC_ref->{$keyBases};
								}
							}
						}
					}
					else {$bases=uc(substr($_[5],$pos,1))}
					if (length($bases)==1) {
						if ($bases ne $refBase) {
							my $lpval=0;
							for (my $i=0;$i<$countDiffBases;$i++) {$lpval+=&TestFisher($countDiffBases-$i,$countRefBases+$i,0+$i,($countDiffBases+$countRefBases)-$i)}
							my $pval=1-$lpval;
							if ($pval<=0) {$pval=&TestFisher($countDiffBases,$countRefBases,0,$countDiffBases+$countRefBases)}
							if ($pval<=$maxPval) {($var_ref,$bases)=&printVCF($_[6],$pos,$refBase,$bases,$countBasesW_ref,$countQual_ref,$var_ref,$pval,$countBasesC_ref)}
							else {$bases=uc(substr($_[5],$pos,1))}
							if ($bases eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
							elsif ($bases eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
							else {}
						}
						else {
							if ($bases eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
							elsif ($bases eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
							else {}
						}
					}
					else {
						my $lpval=0;
						for (my $i=0;$i<$countDiffBases;$i++) {$lpval+=&TestFisher($countDiffBases-$i,$countRefBases+$i,0+$i,($countDiffBases+$countRefBases)-$i)}
						my $pval=1-$lpval;
						if ($pval<=0) {$pval=&TestFisher($countDiffBases,$countRefBases,0,$countDiffBases+$countRefBases)}
						if ($pval<=$maxPval) {($var_ref,$bases)=&printVCF($_[6],$pos,$refBase,$bases,$countBasesW_ref,$countQual_ref,$var_ref,$pval,$countBasesC_ref)}
						else {$bases=uc(substr($_[5],$pos,1))}
						if ($bases=~/C/ and $bases=~/G/) {
							$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT");
							$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA");
						}
						elsif ($bases=~/C/) {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
						elsif ($bases=~/G/) {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
						else {}
					}
				}
				else {
					if ($refBase eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
					elsif ($refBase eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
					else {}
					{lock($uncheckedDepthSNV);
					$uncheckedDepthSNV++;}
				}
			}
			elsif ($_[0]->{$pos} and $minSNVperc<1) {
				if (length($_[0]->{$pos})>=$minDepthSNV) {
					($countBasesW_ref,$countQual_ref,$total)=&countBasesSingle($countBasesW_ref,$countQual_ref,$total,$_[0]->{$pos},$_[3]->{$pos},"W");
					my $bases="";
					my $countDiffBases=0;
					my $countRefBases=0;
					if ($total>0) {
						for my $keyBases (keys %{$countBasesW_ref}) {
							if (($countBasesW_ref->{$keyBases}/$total)>=$minSNVperc) {
								$bases.=$keyBases;
								if ($keyBases ne $refBase) {$countDiffBases+=$countBasesW_ref->{$keyBases}}
								else {$countRefBases+=$countBasesW_ref->{$keyBases}}
							}
						}
					}
					else {$bases=uc(substr($_[5],$pos,1))}
					if (length($bases)==1) {
						if ($bases ne $refBase) {
							my $lpval=0;
							for (my $i=0;$i<$countDiffBases;$i++) {$lpval+=&TestFisher($countDiffBases-$i,$countRefBases+$i,0+$i,($countDiffBases+$countRefBases)-$i)}
							my $pval=1-$lpval;
							if ($pval<=0) {$pval=&TestFisher($countDiffBases,$countRefBases,0,$countDiffBases+$countRefBases)}
							if ($pval<=$maxPval) {($var_ref,$bases)=&printVCF($_[6],$pos,$refBase,$bases,$countBasesW_ref,$countQual_ref,$var_ref,$pval,$countBasesC_ref)}
							else {$bases=uc(substr($_[5],$pos,1))}
							if ($bases eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
							else {}
						}
						else {
							if ($bases eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
							else {}
						}
					}
					else {
						my $lpval=0;
						for (my $i=0;$i<$countDiffBases;$i++) {$lpval+=&TestFisher($countDiffBases-$i,$countRefBases+$i,0+$i,($countDiffBases+$countRefBases)-$i)}
						my $pval=1-$lpval;
						if ($pval<=0) {$pval=&TestFisher($countDiffBases,$countRefBases,0,$countDiffBases+$countRefBases)}
						if ($pval<=$maxPval) {($var_ref,$bases)=&printVCF($_[6],$pos,$refBase,$bases,$countBasesW_ref,$countQual_ref,$var_ref,$pval,$countBasesC_ref)}
						else {$bases=uc(substr($_[5],$pos,1))}
						if ($bases=~/C/) {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
						else {}
					}
				}
				else {
					if ($refBase eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}	
					else {}
					{lock($uncheckedDepthSNV);
					$uncheckedDepthSNV++;}
				}
			}
			elsif ($_[1]->{$pos} and $minSNVperc<1) {
				if (length($_[1]->{$pos})>=$minDepthSNV) {
					($countBasesC_ref,$countQual_ref,$total)=&countBasesSingle($countBasesC_ref,$countQual_ref,$total,$_[1]->{$pos},$_[4]->{$pos},"C");
					my $bases="";
					my $countDiffBases=0;
					my $countRefBases=0;
					if ($total>0) {
						for my $keyBases (keys %{$countBasesC_ref}) {
							if (($countBasesC_ref->{$keyBases}/$total)>=$minSNVperc) {
								$bases.=$keyBases;
								if ($keyBases ne $refBase) {$countDiffBases+=$countBasesC_ref->{$keyBases}}
								else {$countRefBases+=$countBasesC_ref->{$keyBases}}
							}
						}
					}
					else {$bases=uc(substr($_[5],$pos,1))}
					if (length($bases)==1) {
						if ($bases ne $refBase) {
							my $lpval=0;
							for (my $i=0;$i<$countDiffBases;$i++) {$lpval+=&TestFisher($countDiffBases-$i,$countRefBases+$i,0+$i,($countDiffBases+$countRefBases)-$i)}
							my $pval=1-$lpval;
							if ($pval<=0) {$pval=&TestFisher($countDiffBases,$countRefBases,0,$countDiffBases+$countRefBases)}
							if ($pval<=$maxPval) {($var_ref,$bases)=&printVCF($_[6],$pos,$refBase,$bases,$countBasesW_ref,$countQual_ref,$var_ref,$pval,$countBasesC_ref)}
							else {$bases=uc(substr($_[5],$pos,1))}
							if ($bases eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
							else {}
						}
						else {
							if ($bases eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
							else {}
						}
					}
					else {
						my $lpval=0;
						for (my $i=0;$i<$countDiffBases;$i++) {$lpval+=&TestFisher($countDiffBases-$i,$countRefBases+$i,0+$i,($countDiffBases+$countRefBases)-$i)}
						my $pval=1-$lpval;
						if ($pval<=0) {$pval=&TestFisher($countDiffBases,$countRefBases,0,$countDiffBases+$countRefBases)}
						if ($pval<=$maxPval) {($var_ref,$bases)=&printVCF($_[6],$pos,$refBase,$bases,$countBasesW_ref,$countQual_ref,$var_ref,$pval,$countBasesC_ref)}
						else {$bases=uc(substr($_[5],$pos,1))}
						if ($bases=~/G/) {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
						else {}
					}
				}
				else {
					if ($refBase eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
					else {}
					{lock($uncheckedDepthSNV);
					$uncheckedDepthSNV++;}
				}	
			}
			else { # without data from both strands
				if ($refBase eq "C") {$meth_ref=&methylationCall($_[0]->{$pos},$_[3]->{$pos},$meth_ref,$pos,"CT")}
				elsif ($refBase eq "G") {$meth_ref=&methylationCall($_[1]->{$pos},$_[4]->{$pos},$meth_ref,$pos,"GA")}
				else {}
			}
		}
		&printMeth($meth_ref,$var_ref,$_[5],$context,$_[6]);
	}
}

#countBases
#0->hash variation 1->hash quals 2->total count 3->bases in position 4->quality in position 
#5->strand
sub countBases {
	if ($_[5] eq "W") {
		for (my $i=0;$i<=length($_[3]);$i++) {
			if (substr($_[3],$i,1) eq "A") {
				$_[0]->{"A"}+=1;
				$_[1]->{"A"}.=substr($_[4],$i,1);
				$_[2]++;
			}
			if (substr($_[3],$i,1) eq "G") {
				$_[0]->{"G"}+=1;
				$_[1]->{"G"}.=substr($_[4],$i,1);
				$_[2]++;
			}
		}
	}
	if ($_[5] eq "C") {
		for (my $i=0;$i<=length($_[3]);$i++) {
			if (substr($_[3],$i,1) eq "C") {
				$_[0]->{"C"}+=1;
				$_[1]->{"C"}.=substr($_[4],$i,1);
				$_[2]++;
			}
			if (substr($_[3],$i,1) eq "T") {
				$_[0]->{"T"}+=1;
				$_[1]->{"T"}.=substr($_[4],$i,1);
				$_[2]++;
			}
		}
	}
	return($_[0],$_[1],$_[2]);
}

#countBases Correction
#0->hash variationW 1->totalBis 2->total count 3->bases in position 4->quality in position 
#5->strand 6->hash variationC
sub countBasesCorrect {
	if ($_[5] eq "W") {
		if ($_[6]->{"C"}/$_[2]>=$minSNVperc and $_[6]->{"T"}/$_[2]>=$minSNVperc) {
			my $bisBases=0;
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "C") {$bisBases++}
				elsif (substr($_[3],$i,1) eq "T") {$bisBases++}
			}
			$_[0]->{"C"}+=int(($bisBases*($_[6]->{"C"}/$_[2]))+0.5);
			$_[0]->{"T"}+=int(($bisBases*($_[6]->{"T"}/$_[2]))+0.5);
			$_[1]+=int(($bisBases*($_[6]->{"C"}/$_[2]))+0.5)+int(($bisBases*($_[6]->{"T"}/$_[2]))+0.5);
		}
		elsif ($_[6]->{"C"}/$_[2]>=$minSNVperc) {
			my $bisBases=0;
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "C") {$bisBases++}
				elsif (substr($_[3],$i,1) eq "T") {$bisBases++}
			}
			$_[0]->{"C"}+=$bisBases;
			$_[1]+=$bisBases;
		}
		elsif ($_[6]->{"T"}/$_[2]>=$minSNVperc) {
			my $bisBases=0;
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "C") {$bisBases++}
				elsif (substr($_[3],$i,1) eq "T") {$bisBases++}
			}
			$_[0]->{"T"}+=$bisBases;
			$_[1]+=$bisBases;
		}
	}
	if ($_[5] eq "C") {
		if ($_[0]->{"G"}/$_[2]>=$minSNVperc and $_[0]->{"A"}/$_[2]>=$minSNVperc) {
			my $bisBases=0;
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "G") {$bisBases++}
				elsif (substr($_[3],$i,1) eq "A") {$bisBases++}
			}
			$_[6]->{"G"}+=int(($bisBases*($_[0]->{"G"}/$_[2]))+0.5);
			$_[6]->{"A"}+=int(($bisBases*($_[0]->{"A"}/$_[2]))+0.5);
			$_[1]+=$bisBases;
		}
		elsif ($_[0]->{"G"}/$_[2]>=$minSNVperc) {
			my $bisBases=0;
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "G") {$bisBases++}
				elsif (substr($_[3],$i,1) eq "A") {$bisBases++}
			}
			$_[6]->{"G"}+=$bisBases;
			$_[1]+=$bisBases;
		}
		elsif ($_[0]->{"A"}/$_[2]>=$minSNVperc) {
			my $bisBases=0;
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "G") {$bisBases++}
				elsif (substr($_[3],$i,1) eq "A") {$bisBases++}
			}
			$_[6]->{"A"}+=$bisBases;
			$_[1]+=$bisBases;
		}
	}
	return($_[0],$_[6],$_[1]);
}

#countBases single strand
#0->hash variation 1->hash quals 2->total count 3->bases in position 4->quality in position 
#5->strand 6-> reference
sub countBasesSingle {
	if ($_[5] eq "W") {
		my $bisBase=0;
		for (my $i=0;$i<=length($_[3]);$i++) {
			if (substr($_[3],$i,1) eq "C") {$bisBase++}
			if (substr($_[3],$i,1) eq "T") {$bisBase++}
		}
		if ($bisBase/length($_[3])>=$minSNVperc) {}
		else {
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "A") {
					$_[0]->{"A"}+=1;
					$_[1]->{"A"}.=substr($_[4],$i,1);
					$_[2]++;
				}
				if (substr($_[3],$i,1) eq "G") {
					$_[0]->{"G"}+=1;
					$_[1]->{"G"}.=substr($_[4],$i,1);
					$_[2]++;
				}
			}	
		}
	}
	if ($_[5] eq "C") {
		my $bisBase=0;
		for (my $i=0;$i<=length($_[3]);$i++) {
			if (substr($_[3],$i,1) eq "G") {$bisBase++}
			if (substr($_[3],$i,1) eq "A") {$bisBase++}
		}
		if ($bisBase/length($_[3])>=$minSNVperc) {}
		else {
			for (my $i=0;$i<=length($_[3]);$i++) {
				if (substr($_[3],$i,1) eq "C") {
					$_[0]->{"C"}+=1;
					$_[1]->{"C"}.=substr($_[4],$i,1);
					$_[2]++;
				}
				if (substr($_[3],$i,1) eq "T") {
					$_[0]->{"T"}+=1;
					$_[1]->{"T"}.=substr($_[4],$i,1);
					$_[2]++;
				}
			}	
		}
	}
	return($_[0],$_[1],$_[2]);
}

#Calling methylation levels
#0-> bases positions 1-> qual positions 2-> methylation hash
#3-> position 4-> methylation bases
sub methylationCall {
	if ($_[0]) {
		my $methValues=0;
		my $totalValues=0;
		my $asciiValues="";
		for (my $i=0;$i<=length($_[0])-1;$i++) {
			if (substr($_[0],$i,1) eq substr($_[4],0,1)) {
				$methValues+=1;
				$totalValues+=1;
				$asciiValues.=substr($_[1],$i,1);
			}
			elsif (substr($_[0],$i,1) eq substr($_[4],1,1)) {
				$totalValues+=1;
				$asciiValues.=substr($_[1],$i,1);
			}
			else {}
		}
		if ($_[2]->{$_[3]}) {
			my @splitData=split(/\t/,$_[2]->{$_[3]});
			$methValues+=$splitData[0];
			$totalValues+=$splitData[1];
			$asciiValues.=$splitData[2];
			if ($totalValues>0) {$_[2]->{$_[3]}="$methValues\t$totalValues\t$asciiValues"}
		}
		else {
			if ($totalValues>0) {$_[2]->{$_[3]}="$methValues\t$totalValues\t$asciiValues"}	
		}
	}
	return($_[2]);
}

#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome 5->pos 6->init base
sub checkCinCG {
	my $sampleContext=$_[6];
	my $forwData=".\t.\t.";
	my $revData=".\t.\t.";
	if ($_[1]->{$_[5]+1}) {
		if ($_[1]->{$_[5]+1}=~/G/i) {
			$sampleContext.=$bases{$_[1]->{$_[5]+1}};
			$forwData=&methFormat($_[0]->{$_[5]},$forwData);
			if ($_[0]->{$_[5]+1}) {$revData=&methFormat($_[0]->{$_[5]+1},$revData);}
			if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
				open OUT, ">>$outDir/CG_$_[4].output";
				my $pos1=$_[5]+1;
				print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
				close OUT;
				my $pos0=$_[5];
				if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CG",$forwData,$revData)}
				if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CG",$forwData,$revData)}
			}
		}
		else {}
	}
	elsif (length($_[2])>$_[5]+1) {
		if (uc(substr($_[2],$_[5]+1,1)) eq "G") {
			$sampleContext.=uc(substr($_[2],$_[5]+1,1));
			$forwData=&methFormat($_[0]->{$_[5]},$forwData);
			if ($_[0]->{$_[5]+1}) {$revData=&methFormat($_[0]->{$_[5]+1},$revData);}
			if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
				open OUT, ">>$outDir/CG_$_[4].output";
				my $pos1=$_[5]+1;
				print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
				close OUT;
				my $pos0=$_[5];
				if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CG",$forwData,$revData)}
				if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CG",$forwData,$revData)}
			}
		}
	}
	else {}
}

#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome 5->pos 6->first base
sub checkGinCG {
	my $sampleContext=$_[6];
	my $forwData=".\t.\t.";
	my $revData=".\t.\t.";
	if ($_[1]->{$_[5]-1}) {
		if ($_[1]->{$_[5]-1}=~/C/i) {
			if (!$_[0]->{$_[5]-1}) {
				$sampleContext=$bases{$_[1]->{$_[5]-1}}.$sampleContext;
				$revData=&methFormat($_[0]->{$_[5]},$revData);
				if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
					open OUT, ">>$outDir/CG_$_[4].output";
					my $pos1=$_[5];
					print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
					close OUT;
					my $pos0=$_[5]-1;
					if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CG",$forwData,$revData)}
					if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CG",$forwData,$revData)}
				}
			}
		}
		else {}
	}
	elsif ($_[5]-1>=0) {
		if (uc(substr($_[2],$_[5]-1,1)) eq "C") {
			if (!$_[0]->{$_[5]-1}) {
				$sampleContext=uc(substr($_[2],$_[5]-1,1)).$sampleContext;
				$revData=&methFormat($_[0]->{$_[5]},$revData);
				if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
					open OUT, ">>$outDir/CG_$_[4].output";
					my $pos1=$_[5];
					print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
					close OUT;
					my $pos0=$_[5]-1;
					if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CG",$forwData,$revData)}
					if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CG",$forwData,$revData)}
				}
			}
		}
		else {}
	}
	else {}
}

#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome 5->pos 6->init base
sub checkCinCHG {
	my $sampleContext=$_[6];
	my $forwData=".\t.\t.";
	my $revData=".\t.\t.";
	
	#Get second base in context
	if ($_[1]->{$_[5]+1}) {$sampleContext.=$bases{$_[1]->{$_[5]+1}}}
	elsif (uc(substr($_[2],$_[5]+1,1))) {$sampleContext.=uc(substr($_[2],$_[5]+1,1))}
	else{}
	
	if ($_[1]->{$_[5]+2}) {
		if ($_[1]->{$_[5]+2}=~/G/i) {
			$sampleContext.=$bases{$_[1]->{$_[5]+2}};
			if (substr($sampleContext,1,1) ne "G") {$forwData=&methFormat($_[0]->{$_[5]},$forwData)}
			if (substr($sampleContext,1,1) ne "C") {
				if ($_[0]->{$_[5]+2}) {$revData=&methFormat($_[0]->{$_[5]+2},$revData);}	
			}
			if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
				open OUT, ">>$outDir/CHG_$_[4].output";
				my $pos1=$_[5]+1;
				print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
				close OUT;
				my $pos0=$_[5];
				if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CHG",$forwData,$revData)}
				if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CHG",$forwData,$revData)}
			}
		}
		else {}
	}
	elsif (length($_[2])>$_[5]+2) {
		if (uc(substr($_[2],$_[5]+2,1)) eq "G") {
			$sampleContext.=uc(substr($_[2],$_[5]+2,1));
			if (substr($sampleContext,1,1) ne "G") {$forwData=&methFormat($_[0]->{$_[5]},$forwData)}
			if (substr($sampleContext,1,1) ne "C") {
				if ($_[0]->{$_[5]+2}) {$revData=&methFormat($_[0]->{$_[5]+2},$revData);}
			}
			if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
				open OUT, ">>$outDir/CHG_$_[4].output";
				my $pos1=$_[5]+1;
				print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
				close OUT;
				my $pos0=$_[5];
				if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CHG",$forwData,$revData)}
				if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CHG",$forwData,$revData)}
			}
		}
	}
	else {}
}

#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome 5->pos 6->first base
sub checkGinCHG {
	my $sampleContext=$_[6];
	my $forwData=".\t.\t.";
	my $revData=".\t.\t.";
	
	#Get second base in context
	if ($_[1]->{$_[5]-1}) {$sampleContext=$bases{$_[1]->{$_[5]-1}}.$sampleContext}
	elsif (uc(substr($_[2],$_[5]-1,1))) {$sampleContext=uc(substr($_[2],$_[5]-1,1)).$sampleContext}
	else{}
	
	if ($_[1]->{$_[5]-2}) {
		if ($_[1]->{$_[5]-2}=~/C/i) {
			if (!$_[0]->{$_[5]-2}) {
				$sampleContext=$bases{$_[1]->{$_[5]-2}}.$sampleContext;
				if (substr($sampleContext,1,1) ne "C") {
					$revData=&methFormat($_[0]->{$_[5]},$revData);
				}
				if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
					open OUT, ">>$outDir/CHG_$_[4].output";
					my $pos1=$_[5]-1;
					print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
					close OUT;
					my $pos0=$_[5]-2;
					if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CHG",$forwData,$revData)}
					if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CHG",$forwData,$revData)}
				}
			}
		}
		else {}
	}
	elsif ($_[5]-2>=0) {
		if (uc(substr($_[2],$_[5]-2,1)) eq "C") {
			if (!$_[0]->{$_[5]-2}) {
				$sampleContext=uc(substr($_[2],$_[5]-2,1)).$sampleContext;
				if (substr($sampleContext,1,1) ne "C") {
					$revData=&methFormat($_[0]->{$_[5]},$revData);
				}
				if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
					open OUT, ">>$outDir/CHG_$_[4].output";
					my $pos1=$_[5]-1;
					print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
					close OUT;
					my $pos0=$_[5]-2;
					if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CHG",$forwData,$revData)}
					if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CHG",$forwData,$revData)}
				}
			}
		}
		else {}
	}
	else {}
}

#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome 5->pos 6->init base
sub checkCinCHH {
	my $sampleContext=$_[6];
	my $forwData=".\t.\t.";
	my $revData=".\t.\t.";
	
	#Get second base in context
	if ($_[1]->{$_[5]+1}) {
		if ($_[1]->{$_[5]+1}=~/[ATC]/) {$sampleContext.=$bases{$_[1]->{$_[5]+1}}}
		else {}
	}
	elsif (length($_[2])>$_[5]+1) {
		if (uc(substr($_[2],$_[5]+1,1)) ne "G") {$sampleContext.=uc(substr($_[2],$_[5]+1,1))}
	}
	else{}
	
	#Get third base in context
	if ($_[1]->{$_[5]+2}) {
		if ($_[1]->{$_[5]+2}=~/[ATC]/) {$sampleContext.=$bases{$_[1]->{$_[5]+2}}}
		else {}
	}
	elsif (length($_[2])>$_[5]+2) {
		if (uc(substr($_[2],$_[5]+2,1)) ne "G") {$sampleContext.=uc(substr($_[2],$_[5]+2,1))}
	}
	else{}
	
	if (length($sampleContext)==3) {
		$forwData=&methFormat($_[0]->{$_[5]},$forwData);
		if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
			open OUT, ">>$outDir/CHH_$_[4].output";
			my $pos1=$_[5]+1;
			print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
			close OUT;
			my $pos0=$_[5];
			if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CHH",$forwData,$revData)}
			if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CHH",$forwData,$revData)}
		}
	}
}

#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome 5->pos 6->init base
sub checkGinCHH {
	my $sampleContext=$_[6];
	my $forwData=".\t.\t.";
	my $revData=".\t.\t.";
	
	#Get second base in context
	if ($_[1]->{$_[5]-1}) {
		if ($_[1]->{$_[5]-1}=~/[ATG]/) {$sampleContext=$bases{$_[1]->{$_[5]-1}}.$sampleContext}
		else {}
	}
	elsif ($_[5]-1>=0) {
		if (uc(substr($_[2],$_[5]-1,1)) ne "C") {$sampleContext=uc(substr($_[2],$_[5]-1,1)).$sampleContext}
	}
	else{}
	
	#Get third base in context
	if ($_[1]->{$_[5]-2}) {
		if ($_[1]->{$_[5]-2}=~/[ATG]/) {$sampleContext=$bases{$_[1]->{$_[5]-2}}.$sampleContext}
		else {}
	}
	elsif ($_[5]-2>=0) {
		if (uc(substr($_[2],$_[5]-2,1)) ne "C") {$sampleContext=uc(substr($_[2],$_[5]-2,1)).$sampleContext}
	}
	else{}
	
	if (length($sampleContext)==3) {
		$revData=&methFormat($_[0]->{$_[5]},$revData);
		if ($forwData ne ".\t.\t." or $revData ne ".\t.\t.") {
			open OUT, ">>$outDir/CHH_$_[4].output";
			my $pos1=$_[5]-1;
			print OUT "$_[4]\t$pos1\t$sampleContext\t$forwData\t$revData\n";
			close OUT;
			my $pos0=$_[5]-2;
			if ($bedOut eq "Y") {&bedFormat($_[4],$pos0,"CHH",$forwData,$revData)}
			if ($wigOut eq "Y") {&wigFormat($_[4],$pos1,"CHH",$forwData,$revData)}
		}
	}
}

#Printing methylation levels
#0-> methylation pos hash 1-> variation pos hash 2-> reference 3-> context 4->chromosome
sub printMeth {
	for my $pos (sort {$a<=>$b} keys %{$_[0]}) {
		if ($_[3] eq "CG" or $_[3] eq "ALL") {
			if ($_[1]->{$pos}) {
				#Both C&G
				if ($_[1]->{$pos}=~/C/i and $_[1]->{$pos}=~/G/i) {
					#Checking C
					&checkCinCG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
					#Checking G
					&checkGinCG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				#Only C
				elsif ($_[1]->{$pos}=~/C/i) {
					&checkCinCG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				#Only G
				elsif ($_[1]->{$pos}=~/G/i) {
					&checkGinCG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				else {}
			}
			else {
				if (uc(substr($_[2],$pos,1)) eq "C") {
					&checkCinCG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,uc(substr($_[2],$pos,1)));
				}
				elsif (uc(substr($_[2],$pos,1)) eq "G") {
					&checkGinCG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,uc(substr($_[2],$pos,1)));
				}
				else {}
			}
		}
		if ($_[3] eq "CHG" or $_[3] eq "ALL") {
			if ($_[1]->{$pos}) {
				#Both C&G
				if ($_[1]->{$pos}=~/C/i and $_[1]->{$pos}=~/G/i) {
					#Checking C
					&checkCinCHG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
					#Checking G
					&checkGinCHG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				#Only C
				elsif ($_[1]->{$pos}=~/C/i) {
					&checkCinCHG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				#Only G
				elsif ($_[1]->{$pos}=~/G/i) {
					&checkGinCHG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				else {}
			}
			else {
				if (uc(substr($_[2],$pos,1)) eq "C") {
					&checkCinCHG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,uc(substr($_[2],$pos,1)));
				}
				elsif (uc(substr($_[2],$pos,1)) eq "G") {
					&checkGinCHG($_[0],$_[1],$_[2],$_[3],$_[4],$pos,uc(substr($_[2],$pos,1)));
				}
				else {}
			}
		}
		if ($_[3] eq "CHH" or $_[3] eq "ALL") {
			if ($_[1]->{$pos}) {
				#Both C&G
				if ($_[1]->{$pos}=~/C/i and $_[1]->{$pos}=~/G/i) {
					#Checking C
					&checkCinCHH($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
					#Checking G
					&checkGinCHH($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				#Only C
				elsif ($_[1]->{$pos}=~/C/i) {
					&checkCinCHH($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				#Only G
				elsif ($_[1]->{$pos}=~/G/i) {
					&checkGinCHH($_[0],$_[1],$_[2],$_[3],$_[4],$pos,$bases{$_[1]->{$pos}});
				}
				else {}
			}
			else {
				if (uc(substr($_[2],$pos,1)) eq "C") {
					&checkCinCHH($_[0],$_[1],$_[2],$_[3],$_[4],$pos,uc(substr($_[2],$pos,1)));
				}
				elsif (uc(substr($_[2],$pos,1)) eq "G") {
					&checkGinCHH($_[0],$_[1],$_[2],$_[3],$_[4],$pos,uc(substr($_[2],$pos,1)));
				}
				else {}
			}
		}
	}
}

#Format Output
#0-> meth values 1-> output string
sub methFormat {
	my @split=split(/\t/,$_[0]);
	if ($split[1]>=$minDepthMeth) {
		my $outQ=0;
		for (my $i=0;$i<=length($split[2])-1;$i++) {$outQ+=&{$subQscore{$qscore}}(substr($split[2],$i,1))}
		$outQ=int(($outQ/length($split[2]))+0.5);
		$_[1]="$split[0]\t$split[1]\t$outQ";
	}
	else {
		{lock($discardDepthPos);
		$discardDepthPos++;}
	}
	return ($_[1]);
}

#BED format output
#0-> chromosome 1-> position 2-> context 3-> watson methylation 4-> crick methylation
sub bedFormat {
	my @splitW=split(/\t/,$_[3]);
	my @splitC=split(/\t/,$_[4]);
	my $chromEnd=$_[1]+length($_[2]);
	my $methValues=0;
	my $coverage=0;
	if ($splitW[0] ne ".") {
		$methValues+=$splitW[0];
		$coverage+=$splitW[1];
	}
	if ($splitC[0] ne ".") {
		$methValues+=$splitC[0];
		$coverage+=$splitC[1];
	}
	my $ratio=int((($methValues/$coverage)*1000)+0.5);
	open BEDOUT, ">>$outDir/$_[2]_$_[0].bed";
	print BEDOUT "$_[0]\t$_[1]\t$chromEnd\t.\t$ratio\t+\n";
	close BEDOUT;
}

#WIG format output
#0-> chromosome 1-> position 2-> context 3-> watson methylation 4-> crick methylation
sub wigFormat {
	my @splitW=split(/\t/,$_[3]);
	my @splitC=split(/\t/,$_[4]);
	my $methValues=0;
	my $coverage=0;
	if ($splitW[0] ne ".") {
		$methValues+=$splitW[0];
		$coverage+=$splitW[1];
	}
	if ($splitC[0] ne ".") {
		$methValues+=$splitC[0];
		$coverage+=$splitC[1];
	}
	my $ratio=int((($methValues/$coverage)*100)+0.5);
	open WIGOUT, ">>$outDir/$_[2]_$_[0].wig";
	print WIGOUT "$_[1]\t$ratio\n";
	close WIGOUT;
}

#Getting LOG FILE information
#0-> line 1-> context 2-> covered context hash 3-> countMeth hash 4-> countUnmeth 5-> total sum hash
sub countValues {
	my @split=split(/\t/,$_[0]);
	my $meth=0;
	my $cover=0;
	if ($split[3] ne ".") {
		$meth+=$split[3];
		$cover+=$split[4];
	}
	if ($split[6] ne ".") {
		$meth+=$split[6];
		$cover+=$split[7];
	}
	$_[2]->{$split[0]}+=1;
	$_[5]->{$split[0]}{$_[1]}+=($meth/$cover);
	if (($meth/$cover)<=0.2) {$_[4]->{$split[0]}{$_[1]}+=1}
	elsif (($meth/$cover)>=0.8) {$_[3]->{$split[0]}{$_[1]}+=1}
	else {}
	return($_[2],$_[3],$_[4],$_[5]);
}
#0->context 1->coverContext 2->numContext 3->countMeth_ref 4->countUnmeth_ref 5->sumTotal_ref
sub outMethResults {
	for my $keyChrom (keys %{$_[1]}) {
		print LOG " $_[0] methylation on $keyChrom:\n";
		my $name="$keyChrom"."_$_[0]";
		print LOG "  Reference contexts: $_[2]->{$name}\n";
		if (!$$_[2]{$name}) {
			my $cover=sprintf(sprintf("%.2f", ($_[1]->{$keyChrom}/$_[2]->{$name})*100));
			print LOG "  Methylation context coverage: $cover\n";
		}
		else {print LOG "  Methylation context coverage: NA\n";}
		if (!$$_[1]{$keyChrom}) {
			my $mean=sprintf(sprintf("%.2f", ($_[5]->{$keyChrom}{$_[0]}/$_[1]->{$keyChrom})));
			print LOG "  Mean methylation: $mean\n";
			my $percMeth=0;
			if ($_[3]->{$keyChrom}{$_[0]}) {
				$percMeth=sprintf(sprintf("%.2f", ($_[3]->{$keyChrom}{$_[0]}/$_[1]->{$keyChrom})*100));
			}
			print LOG "  Percentage of methylated context: $percMeth\n";
			my $percUnmeth=0;
			if ($_[4]->{$keyChrom}{$_[0]}) {
				$percUnmeth=sprintf(sprintf("%.2f", ($_[4]->{$keyChrom}{$_[0]}/$_[1]->{$keyChrom})*100));
			}
			print LOG "  Percentage of unmethylated context: $percUnmeth\n";
		}
		else {
			print LOG "  Mean methylation: NA\n";
			print LOG "  Percentage of methylated context: NA\n";
			print LOG "  Percentage of unmethylated context: NA\n";
		}
	}
}


sub TestFisher{
	my $overA = &BinomialCoeff(($_[0]+$_[1]),$_[0]);
	my $overB= 	&BinomialCoeff(($_[2]+$_[3]),$_[2]);
	my $under= 	&BinomialCoeff(($_[0]+$_[1]+$_[2]+$_[3]),($_[0]+$_[2]));
	#my $Coeff= ($overA*$overB) / $under;
	my $Coeff= exp($overA+$overB-$under);
	return ($Coeff);
}

## 0->over; 1-->under
#sub BinomialCoeff {
#	my $r=1;
#	for (1 .. $_[1]) {	$r *= $_[0] + 1 - $_, $r /= $_ }
#	return($r)
#}

## 0->over; 1-->under
sub BinomialCoeff	{
	my $over = $_[0];
	my $under = $_[1];
	if($under > $over){
		return -1;
	}
	elsif($under == 0){
		return 0;
	}
	elsif($under < $over - $under){
		# first sum
		my $first = 0;
		my $second = 0;
		for(my $i = $over-$under+1; $i <= $over; $i++){
			$first += log($i);
		}
		# $second sum
		for(my $i = 1; $i <= $under; $i++){
			$second += log($i);
		}
		return $first-$second;
	}
	else{
		my $first = 0;
		my $second = 0;
		for(my $i = $under+1; $i <= $over; $i++){
			$first +=  log($i);
		}
		# second sum
		for(my $i = 1; $i <= $over - $under; $i++){
			$second +=  log($i);
		}
		return $first-$second;	
	}
}
