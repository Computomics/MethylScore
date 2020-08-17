#!/usr/bin/env perl
# Copyright (C) 2016-2018 Computomics GmbH

use Getopt::Long;
use Thread::Pool;
use Thread::Conveyor::Array;
use Thread::Conveyor::Throttled;
use Thread::Conveyor::Tied;
use Thread::Tie::Array;
use Config::Simple;
use FindBin qw($Bin);
use Cwd 'abs_path';
use File::Path qw(remove_tree make_path);
use File::Copy;
use File::Tee 'tee';
use File::Which;
use File::Basename;
use POSIX qw(strftime);
use List::Util qw(min);
use Data::Dumper;

use strict;
use experimental qw(smartmatch);



##########
# DISCLAIMER: Pipeline starts for now after the mapping stage of
#             each individual sample (including replicates), thus
#             some samples can have multiple reps, some none


# Notes: o  all CONFIG variables are in scope $CFG:: and uppercase,
#        o  all user input start with $OPT_


###########################################################
### Global variables:
###########################################################

our $version = "0.1.17";

our %STAGES = (
  dedup =>      {
                  order => 1,
                  parallel => 1,
                  time => "10:00:00",
                  mem => "4G",
                  active => 1,
                  dependson => undef,
                  progressbar => 1,
                  needs_mappingfiles => 1,
                  descr => "Processing mapping files"
                },
  # read_stats => {
  #                  order => 2,
  #                  parallel => 1,
  #                  active => 0,
  #                  dependson => \$CFG::STATISTICS
  #               },
  # cov_stats =>  {
  #                  order => 3,
  #                  parallel => 1,
  #                  active => 0,
  #                  dependson => \$CFG::STATISTICS
  #               },
  consensus =>  {
                  order => 4,
                  parallel => 1,
                  time => "40:00:00",
                  mem => "4G",
                  active => 1,
                  dependson => undef,
                  progressbar => 1,
                  needs_mappingfiles => 1,
                  descr => "Consensus calling"
                },
  matrix =>     {
                  order => 5,
                  parallel => 1,
                  time => "10:00:00",
                  mem => "4G",
                  active => 1,
                  dependson => undef,
                  progressbar => 1,
                  needs_mappingfiles => 0,
                  descr => "Generating chromosomal genome matrices"
                },
  matrixWG =>  {
                 order => 6,
                 parallel => 0,
                 time => "10:00:00",
                 mem => "4G",
                 active => 1,
                 dependson => undef,
                 progressbar => 0,
                 needs_mappingfiles => 0,
                 descr => "Generating global genome matrix"
               },
  MRs =>       {
                 order => 7,
                 parallel => 1,
                 time => "10:00:00",
                 mem => "4G",
                 active => 1,
                 dependson => undef,
                 progressbar => 1,
                 needs_mappingfiles => 0,
                 descr => "Determining methylated regions"
               },
  igv =>       {
                 order => 8,
                 parallel => 0,
                 time => "10:00:00",
                 mem => "4G",
                 active => 0,
                 dependson => undef,
                 progressbar => 0,
                 needs_mappingfiles => 0,
                 descr => "Generating IGV file"
               },
  DMRs =>      {
                 order => 9,
                 parallel => 1,
                 time => "40:00:00",
                 mem => "4G",
                 active => 1,
                 dependson => undef,
                 progressbar => 1,
                 needs_mappingfiles => 0,
                 descr => "Determining differentially methylated regions"
               },
  # region_stats => {
  #                    order => 10,
  #                    parallel => 0,
  #                    active => 0,
  #                    dependson => \$CFG::STATISTICS
  #                 }
);

# user-specified parameters:
our $OPT_config_file = "";
our $OPT_stages;# = join(",", keys %STAGES);
our $OPT_project_folder;
our $OPT_sample_sheet;
our %OPT_contexts = ();
our $OPT_threads;
our $OPT_run_local = 0;
our $OPT_force_rerun = 0;
our $OPT_rm_intmd;
our $OPT_donot_dedup;
our $OPT_verbose = 1;
our $SGE = 0;
our $DEBUG = 0;

# global data structures:
our %SAMPLES = ();
# key: sample ID
# --- $ idx: index of the sample
# --- @ files: contains file names of all technical replicates
#              (at the moment: mapping files)
# --- $ reference: file name of sample-specific reference genome fasta file
# --- $ chrs: number of chromosomes for which there are mappings in bam file

our @CHROM_ORDER = ();
# list of chromosome names in sample reference files (same number and order of chroms
# in all samples required!)

our %JOBS_DEDUP = ();
# key: JobIDs of all queued or running JOBS_DEDUP (finished ones are deleted)
# value: sampleID

our %JOBS_CONS = ();
# key: JobIDs of all queued or running JOBS_CONS (finished ones are deleted)
# value: chr

our %JOBS_MATRIX = ();
# key: JobIDs of all queued or running JOBS_MATRIX (finished ones are deleted)
# value: chr

our %JOBS_MRs = ();
# key: JobIDs of all queued or running JOBS_MRs (finished ones are deleted)
# value: sampleID

our %CHRS_PER_SAMPLE = ();
# key: sampleID
# value: chromosomes for which there are mappings in .bam file

our %SAMPLES_PER_CHR = ();
# key: chromosomes for which there are mappings in .bam file
# value: number of samples for which consensus is running on that chr

our %UNION_CHROMS = ();
# key: chromosome
# value: array of samples that have consensus on that chromosome

our %MRS_PER_SAMPLE = ();
# key: sampleID
# value: file name of sample-specific MR

our %JOBS_PENDING;
# key: continuous internal Job ID ($JOBID)
# value: %command
#   $ str: command str (to run locally)
#   $ id: same variable as key
#   $ jobid: Job ID from the cluster or from the Thread::Pool

our %JOBS_RUNNING;

our $JOBID=0;

our $ACCOUNTING_FILE;
our %JOBS_FINISHED;
our $ACCOUNTING_FILE_TAIL_LEN = 100;
our $MAX_ACCFILE_TAIL_LEN = 500000;
our $ACCOUNTING_FILE_SLEEP = 2;

our $BAMTOOLS_SPLIT;

our $NEED_MAPPINGFILES;

# for DEBUGGING:
our $SLEEP_TOTAL=0;

# hard-coded:
our $ERRMSG = "\n[MethylScore] ERROR !!! ";
our $WARNMSG = "\n[MethylScore] WARNING ! ";
###########################################################
###########################################################



# get command line options:
my $command = GetCom();

### Read variables from config file:
#my $CFG = new Config::Simple($OPT_config_file)->vars();
die $ERRMSG . "Config file '$OPT_config_file' not found!" if (! -e $OPT_config_file);
our $CFG;
Config::Simple->new($OPT_config_file)->import_names('CFG');
check_config_variables();


parse_stages();


# Initialize project's global log file. Following command prints
# everything that is printed to STDERR to a file as well:
move("$CFG::PROJECT_FOLDER/log", "$CFG::PROJECT_FOLDER/log.". (strftime "%F.%R", localtime));
tee STDERR, ">>$CFG::PROJECT_FOLDER/log";

### Log start
if ($OPT_verbose) {
  print STDERR "\n" . "##### Starting MethylScore ##### (" . localtime() . ")\n\nexecuting: $command\n\n";

  print STDERR "########### PARAMETERS:\n";
  foreach my $varname (sort keys %CFG::) {
    my $param = ${$CFG::{$varname}};
    if (ref($param) eq 'ARRAY') {
      print STDERR "$varname: " . join(",", @$param) . "\n";
    } else {
      print STDERR "$varname: $param\n";
    }
  }
  print STDERR "#######################\n\n";
}

# warnings file:
our $WARNINGS = 0;
open WARNINGS, ">$CFG::PROJECT_FOLDER/warnings.txt";

# debug file:
open DEBUG, ">$CFG::PROJECT_FOLDER/debug.out" if ($DEBUG);


###########################################################
# Global variables depending on config file:
our $MAPPING_FOLDER = $CFG::PROJECT_FOLDER . "/01mappings/";
our $CONSENSUS_FOLDER = $CFG::PROJECT_FOLDER . "/02consensus/";
our $MATRIX_FOLDER = $CFG::PROJECT_FOLDER . "/03matrix/";
our $MR_FOLDER = $CFG::PROJECT_FOLDER . "/04MRs/";
our $DMR_FOLDER = $CFG::PROJECT_FOLDER . "/05DMRs/";
our $IGV_FOLDER = $CFG::PROJECT_FOLDER . "/igv";
our $TMP_FOLDER = $CFG::PROJECT_FOLDER . "/tmp";

# Defining folder structure:
mkdir($MAPPING_FOLDER) if (! -e $MAPPING_FOLDER);
mkdir($CONSENSUS_FOLDER) if (! -e $CONSENSUS_FOLDER);
mkdir($MATRIX_FOLDER) if (! -e $MATRIX_FOLDER);
our $matrixfile = $MATRIX_FOLDER."genome_matrix.tsv";
mkdir($MR_FOLDER) if (! -e $MR_FOLDER);
mkdir($DMR_FOLDER) if (! -e $DMR_FOLDER);
mkdir($TMP_FOLDER) if (! -e $TMP_FOLDER);
if ($CFG::CLUSTER && $SGE) {
  # grab cluster accounting file using environment variables:
  die $ERRMSG . "Environment variable SGE_ROOT not set! Cannot find cluster accounting file\n" if (!defined $ENV{"SGE_ROOT"});
  die $ERRMSG . "Environment variable SGE_CELL not set! Cannot find cluster accounting file\n" if (!defined $ENV{"SGE_CELL"});
  $ACCOUNTING_FILE = $ENV{"SGE_ROOT"} . "/" . $ENV{"SGE_CELL"} . "/common/accounting";
  die $ERRMSG . "Cannot find cluster accounting file $ACCOUNTING_FILE\n" if (! -e $ACCOUNTING_FILE);
}

# Initialize threading container that regulates max. no.
# of running threads automatically:
our $Thread_pool;
if (!$CFG::CLUSTER) {
  $Thread_pool = Thread::Pool->new(
    {
      #do => sub { print "--- @_\n"; return 1; },
      do => sub {
print DEBUG "now pooling: " . join(" ", @_) . "\n" if ($DEBUG);
        return system("@_");
      },
      workers => $CFG::THREADS,
      maxjobs => undef,
      minjobs => 1,
    }
  );
}
###########################################################



# read in sample file with mapping files:
read_sample_file();


pipeline();


print STDERR "There are warnings. Check $CFG::PROJECT_FOLDER/warnings.txt!\n" if ($WARNINGS);
exit 0;

# STEP -- 1 -- Merge and deduplicate
#
#   -----------------  -----------------           -----------------
#   | sample1, rep1 |  | sample1, rep2 |           |    sample2    |
#   -----------------  -----------------           -----------------
#             v     MERGE     v                            |
#             -----------------                            |
#             |    sample1    |                            |
#             -----------------                            |
#                     v            DEDUPLICATE             v
#             -----------------                   -----------------
#             |    sample1    |                   |    sample2    |
#             -----------------                   -----------------

# STEP -- 2 -- Run consensus
#
#             -----------------                   -----------------
#             |    sample1    |                   |    sample2    |
#             -----------------                   -----------------
#                     v           split into chrs         v
#             ---- ---- ---- --                   ---- ---- ---- --
#             |  | |  | |  | ||                   |  | |  | |  | ||
#             ---- ---- ---- --                   ---- ---- ---- --
#                     v           methylextract           v
#             ---- ---- ---- --                   ---- ---- ---- --
#             |  | |  | |  | ||                   |  | |  | |  | ||
#             ---- ---- ---- --                   ---- ---- ---- --

# STEP -- 3 -- Generate genome matrix for each chromosome
#
#                consensus                            consensus
#             ---- ---- ---- --                   ---- ---- ---- --
#             |  | |  | |  | ||                   |  | |  | |  | ||
#             ---- ---- ---- --                   ---- ---- ---- --
#                     v       generate_genome_matrix      v
#             ---- ---- ---- --                   ---- ---- ---- --
#             |  | |  | |  | ||                   |  | |  | |  | ||
#             ---- ---- ---- --                   ---- ---- ---- --

# STEP -- 4 -- Call methylated regions (MRs)
#                  sample1                            sample2
#             ---- ---- ---- --                  ---- ---- ---- --
#             |  | |  | |  | || matrices per chr |  | |  | |  | ||
#             ---- ---- ---- --                  ---- ---- ---- --
#                              \   sort into    /
#                               \              /
#                          multi-sample matrix per chr
#                               ---- ---- ---- --
#                               |  | |  | |  | ||
#                               ---- ---- ---- --
#                                       v           merge
#                               -----------------
#                               | genome matrix |
#                               -----------------
#                                v    v    v        call MRs per sample
#                               ---- ---- ----
#                               |S1| |S2| |S3|...   MRs
#                               ---- ---- ----
#
# STEP -- 5 -- Call differentially methylated regions (DMRs)


###########################################################
### Check/Update some variables:
sub check_config_variables {

  abort_pipeline("Cannot find config file!") unless -e $OPT_config_file;

  # folder of config file:
  my $config_folder = abs_path(dirname($OPT_config_file));

  # folder of current working dir:
  my $cwd = Cwd::cwd();

  # Compare following paths from config file with those from this program's argument
  # and in case they're different, take the program's argument one
  my @var_names = ('PROJECT_FOLDER', 'SAMPLE_SHEET', 'BIN_PATH', 'EXTBIN_PATH', 'SCRIPT_PATH',
    'ROI', 'MR_PARAMS');
  my @opt_vars = ($OPT_project_folder, $OPT_sample_sheet, $CFG::BIN_PATH, $CFG::EXTBIN_PATH,
    $CFG::SCRIPT_PATH, $CFG::ROI, $CFG::MR_PARAMS);
  my @cfg_vars = (\$CFG::PROJECT_FOLDER, \$CFG::SAMPLE_SHEET, \$CFG::BIN_PATH,
    \$CFG::EXTBIN_PATH, \$CFG::SCRIPT_PATH, \$CFG::ROI, \$CFG::MR_PARAMS);

#  if ($CFG::BAMTOOLS =~ /\//) {
#    push @var_names, "BAMTOOLS";
#    push @opt_vars, $CFG::BAMTOOLS;
#    push @cfg_vars, \$CFG::BAMTOOLS;
#  }

  my $overwritten_params = "\nSome parameters in config file differ from user input. Will use user input for: ";
  for (my $i=0; $i<@opt_vars; ++$i) {

    # 1. convert all paths into absolute paths:

    #   if relative path on command line given, assume base folder is cwd (for input data),
    #   or relative to MethylScore executable (for binary paths):
    if ($opt_vars[$i] ne "" && $opt_vars[$i] !~ /^\//) {
      if ($var_names[$i] !~ /_PATH$/) {
        $opt_vars[$i] = abs_path($cwd . "/" . $opt_vars[$i]);
      }
      else {
        $opt_vars[$i] = abs_path($Bin . "/../" . $opt_vars[$i]);
      }
    }
    #   if relative path on command line given, assume base folder is config folder (input data),
    #   or relative to MethylScore executable (for binary paths):
    if (${$cfg_vars[$i]} ne "" && ${$cfg_vars[$i]} !~ /^\//) {
      if ($var_names[$i] !~ /_PATH$/) {
        ${$cfg_vars[$i]} = abs_path($config_folder . "/" . ${$cfg_vars[$i]});
      }
      else {
        ${$cfg_vars[$i]} = abs_path($Bin . "/../" . ${$cfg_vars[$i]});
      }
    }

    # 2. override cfg_vars with opt_vars if applicable:
    if ($opt_vars[$i] ne "" && $opt_vars[$i] ne ${$cfg_vars[$i]}) {
      $overwritten_params .= $var_names[$i] . "=" . $opt_vars[$i] . " ";

      ${$cfg_vars[$i]} = $opt_vars[$i];
    }

    # 3. check if folders exist or can be created:
    if (${$cfg_vars[$i]} ne "" && ! -e ${$cfg_vars[$i]}) {
      if ($i == 0) { # PROJECT_FOLDER should be created if not existant
        mkdir(${$cfg_vars[$i]})
          or abort_pipeline("Cannot create ".$var_names[$i].": ".${$cfg_vars[$i]});
      }
      else {
        # the others should exist:
        abort_pipeline("Cannot find $var_names[$i], path does not exist: ".
        ${$cfg_vars[$i]}."\n\nMake sure to provide non-absolute paths relative to the folder that contains the config file");
      }
    }

    # 4. output overwritten params (further below)

  }

  # compare further options between user input and config file:
  if ($CFG::FORCE_RERUN != $OPT_force_rerun) {
    $overwritten_params .= "FORCE_RERUN=$OPT_force_rerun ";
    $CFG::FORCE_RERUN = $OPT_force_rerun;
  }
  if (defined $OPT_threads && $CFG::THREADS != $OPT_threads) {
    $overwritten_params .= "THREADS=$OPT_threads ";
    $CFG::THREADS = $OPT_threads;
  }
  if (defined $OPT_rm_intmd && $CFG::REMOVE_INTMED_FILES != $OPT_rm_intmd) {
    $overwritten_params .= "REMOVE_INTMED_FILES=$OPT_rm_intmd ";
    $CFG::REMOVE_INTMED_FILES = $OPT_rm_intmd;
  }
  if (defined $OPT_donot_dedup && $CFG::DO_DEDUP == $OPT_donot_dedup) {
    $CFG::DO_DEDUP = 1 - $OPT_donot_dedup;
    $overwritten_params .= "DO_DEDUP=$CFG::DO_DEDUP ";
  }

  # output parameters overwritten by user input (4. from above):
  print STDERR $overwritten_params."\n" if ($overwritten_params =~ m/=/);



  # adapt cluster variables to user input:
  if ($OPT_run_local) {
    $CFG::CLUSTER = 0;
  }


  # check BAMTOOLS and SAMTOOLS
  my @vars =  (\$CFG::BAMTOOLS, \$CFG::SAMTOOLS);
  my @execs = ("bamtools", "samtools");
  for (my $i=0; $i<scalar(@vars); ++$i) {
    my $ret = system(${$vars[ $i ]}."> /dev/null 2> /dev/null");
    if ($ret == 1) {
      $ret = system($CFG::EXTBIN_PATH."/".$execs[$i]."> /dev/null 2> /dev/null");
      if ($ret == 1) {
        ${$vars[$i]} = which(${$vars[$i]});
        if (!defined ${$vars[$i]}) {
          abort_pipeline("Cannot find or execute ".$execs[$i].". Please place it in folder ".
            "$CFG::EXTBIN_PATH or install it globally\n");
        }
      }
      else {
        ${$vars[$i]} = $CFG::EXTBIN_PATH."/".$execs[$i];
      }
    }
  }

  $BAMTOOLS_SPLIT = $CFG::BAMTOOLS . ",split";


  # Check DMR parameters:
  die $ERRMSG . "MR_FREQ_CHANGE must be between 0 and 100\n" if ($CFG::MR_FREQ_CHANGE !~ /^\d+$/ || $CFG::MR_FREQ_CHANGE > 100);
  die $ERRMSG . "MR_FREQ_DISTANCE must be a positive number\n" if ($CFG::MR_FREQ_DISTANCE !~
    /^\d+$/);
  die $ERRMSG . "SLIDING_WINDOW_SIZE must be a positive number\n" if ($CFG::SLIDING_WINDOW_SIZE
    !~ /^\d+$/);
  die $ERRMSG . "SLIDING_WINDOW_STEP must be a positive number\n" if ($CFG::SLIDING_WINDOW_STEP
    !~ /^\d+$/);
  die $ERRMSG . "DMR_MIN_C must be a positive number\n" if ($CFG::DMR_MIN_C !~ /^\d+$/);
  die $ERRMSG . "DMR_MIN_COV must be a positive number\n" if ($CFG::DMR_MIN_COV !~ /^\d+$/);
  die $ERRMSG . "MR_BATCH_SIZE must be a positive number\n" if ($CFG::MR_BATCH_SIZE !~ /^\d+$/);
  die $ERRMSG . "HDMR_FOLD_CHANGE must be a positive number\n" if ($CFG::HDMR_FOLD_CHANGE !~
    /^\d+$/);
  die $ERRMSG . "FDR_CUTOFF must be between 0 and 1\n" if ($CFG::FDR_CUTOFF !~ /^\d+\.?\d*$/);

  # check Consensus parameters:
  die $ERRMSG . "MIN_QUAL must be a positive number (max: 40)\n" if ($CFG::MIN_QUAL !~ /^\d+$/ || $CFG::MIN_QUAL > 40);
  die $ERRMSG . "IGNORE_LAST_BP must be a positive number\n" if ($CFG::IGNORE_LAST_BP !~ /^\d+$/);
  die $ERRMSG . "IGNORE_FIRST_BP must be a positive number\n" if ($CFG::IGNORE_FIRST_BP !~ /^\d+$/);

  # check MR parameters:
  die $ERRMSG . "MIN_COVERAGE must be a positive number\n" if ($CFG::MIN_COVERAGE !~ /^\d+$/);
  die $ERRMSG . "DESERT_SIZE must be a positive number\n" if ($CFG::DESERT_SIZE !~ /^\d+$/);
  die $ERRMSG . "MERGE_DIST must be a positive number\n" if ($CFG::MERGE_DIST !~ /^\d+$/);
  die $ERRMSG . "MR_MIN_C must be a positive number or -1\n" if ($CFG::MR_MIN_C !~ /-1|^\d+$/);
  die $ERRMSG . "TRIM_METHRATE must be between 0 and 100\n" if ($CFG::TRIM_METHRATE !~ /^\d+$/ ||
    $CFG::TRIM_METHRATE > 100);

  # check contexts to analyze
  if ($CFG::DMRS_PER_CONTEXT) {
    my @ar = ();
    if (ref($CFG::DMR_CONTEXTS) ne 'ARRAY') {
      push @ar, $CFG::DMR_CONTEXTS;
    }
    else {
      @ar = @{$CFG::DMR_CONTEXTS};
    }

    foreach my $c (@ar) {
      $OPT_contexts{$c} = 1;
      if ($c eq "CG") {
        die $ERRMSG . "CLUSTER_MIN_METH_DIFF_$c must be between 0 and 100\n" if
          ($CFG::CLUSTER_MIN_METH_DIFF_CG !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_DIFF_CG > 100);
        die $ERRMSG . "CLUSTER_MIN_METH_$c must be between 0 and 100\n" if
          ($CFG::CLUSTER_MIN_METH_CG !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_CG > 100);
      } elsif ($c eq "CHG") {
        die $ERRMSG . "CLUSTER_MIN_METH_DIFF_$c must be between 0 and 100\n" if
          ($CFG::CLUSTER_MIN_METH_DIFF_CHG !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_DIFF_CHG > 100);
        die $ERRMSG . "CLUSTER_MIN_METH_$c must be between 0 and 100\n" if
          ($CFG::CLUSTER_MIN_METH_CHG !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_CHG > 100);
      } else {
        die $ERRMSG . "CLUSTER_MIN_METH_DIFF_$c must be between 0 and 100\n" if
          ($CFG::CLUSTER_MIN_METH_DIFF_CHH !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_DIFF_CHH > 100);
        die $ERRMSG . "CLUSTER_MIN_METH_$c must be between 0 and 100\n" if
          ($CFG::CLUSTER_MIN_METH_CHH !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_CHH > 100);
      }
    }
  } else {
    $OPT_contexts{"combined"} = 1;
    die $ERRMSG . "CLUSTER_MIN_METH_DIFF must be between 0 and 100\n" if ($CFG::CLUSTER_MIN_METH_DIFF !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH_DIFF > 100);
    die $ERRMSG . "CLUSTER_MIN_METH must be between 0 and 100\n" if ($CFG::CLUSTER_MIN_METH !~ /^\d+$/ || $CFG::CLUSTER_MIN_METH > 100);
  }


  # take over the cluster resources for each stage:
  my @mem_vars = (\$CFG::DEDUP_MEM, \$CFG::CONSENSUS_MEM, \$CFG::MATRIX_MEM, \$CFG::MATRIXWG_MEM, \$CFG::MRS_MEM, \$CFG::IGV_MEM, \$CFG::DMRS_MEM);
  my @time_vars = (\$CFG::DEDUP_TIME, \$CFG::CONSENSUS_TIME, \$CFG::MATRIX_TIME, \$CFG::MATRIXWG_TIME, \$CFG::MRS_TIME, \$CFG::IGV_TIME, \$CFG::DMRS_TIME);
  my @vars = ("dedup", "consensus", "matrix", "matrixWG", "MRs", "igv", "DMRs");

  for (my $i=0; $i<scalar(@vars); ++$i) {
    abort_pipeline("Format of option " . uc($vars[$i]) . "_MEM not recognized. " .
      "Expected number with trailing M (megabytes) or G (gigabyte)")
      if (${$mem_vars[$i]} !~ /^\d+[GM]$/);
    $STAGES{$vars[$i]}{"mem"} = ${$mem_vars[$i]};

    die $ERRMSG . "Option " . uc($vars[$i]) . "_TIME must be a number!\n"
      if (${$time_vars[$i]} !~ /\d:\d\d:\d\d/);
    $STAGES{$vars[$i]}{"time"} = ${$time_vars[$i]};
  }

}


###########################################################
### Read in samples and their mapping files and write into
### global variable %samples:
sub read_sample_file {

  #TODO MANUAL: no commas in sampleIDs or mapping file names

  # folder of config file:
  my $sheet_folder = abs_path(dirname($CFG::SAMPLE_SHEET));

  print STDERR "[MethylScore] " . localtime() . " - Reading in sample sheet\n" if ($OPT_verbose);
  open SAMPLES, "<", $CFG::SAMPLE_SHEET or abort_pipeline("Cannot open $CFG::SAMPLE_SHEET");

  my $idx = 0;
  my $prev_ref = "";
  while (<SAMPLES>) {

    chomp;
    next if ($_ =~ /^#/);
    my @a = split /,| +|\t/;
    my $sampleID = $a[0];

    if (! exists $SAMPLES{$sampleID}) {
      $SAMPLES{$sampleID}{"idx"} = $idx;
      $idx++;

      my @tmp = ();
      $SAMPLES{$sampleID}{"files"} = \@tmp;
    }

    # NEED_MAPPINGFILES is set in parse_stages()
    if ($NEED_MAPPINGFILES) {
      # make absolute path:
      my $mapping_file = $a[2];
      if ($a[2] !~ /^\//) {
        $mapping_file = abs_path($sheet_folder . "/" . $a[2]);
      }

      abort_pipeline("Cannot find mapping file: " . $a[2]) if (!-e $mapping_file);
      abort_pipeline("Mapping file $mapping_file has no .bam extension!") if ($mapping_file !~ /.bam$/);
      push @{$SAMPLES{$sampleID}{"files"}}, $mapping_file;

      if ($a[1] ne "SE" && $a[1] ne "PE") {
        abort_pipeline("Library type can only be 'SE' or 'PE' (check sample $sampleID)");
      }
      $SAMPLES{$sampleID}{"libtype"} = $a[1];
    }

    my $ref = $a[3];
    if (! defined $ref && $prev_ref ne "") {
      $ref = $prev_ref;
    } else {
      if ($ref !~ /^\//) {
        $ref = abs_path($sheet_folder . "/" . $a[3]);
      }
      else {
        $ref = $a[3];
      }
    }
    abort_pipeline("No reference specified for sample $sampleID") if ($ref eq "");
    abort_pipeline("Cannot find reference file for sample $sampleID") if (! -e $ref);
    if (exists $SAMPLES{$sampleID}{"reference"} && $SAMPLES{$sampleID}{"reference"} ne $ref) {
      abort_pipeline("Different reference files for same sample not supported (sample $sampleID)");
    }
    $SAMPLES{$sampleID}{"reference"} = $ref;
    $prev_ref = $ref;

    my $chroms = `grep "^>" $ref | awk -F ' ' '{print substr(\$1, 2)}'`;
    my @chrom_order = split"\n", $chroms;

    if (scalar(@CHROM_ORDER) == 0 || @chrom_order ~~ @CHROM_ORDER) {}
    else { abort_pipeline("Chromosome order for sample $sampleID differs from the previous
    samples"); }
    @CHROM_ORDER = @chrom_order;

  }
  close SAMPLES;

}


### Set stages inactive that were not selected by user, or that do not meet the crtierion to run (determined by the variable 'dependson')
sub parse_stages {

  $NEED_MAPPINGFILES = 0;
  if ($OPT_stages eq "") {
    $NEED_MAPPINGFILES = 1;
  }
  elsif ($OPT_stages eq "all") {
    foreach (keys %STAGES) {
      $STAGES{$_}{active} = 1;
    }
    $NEED_MAPPINGFILES = 1;
  }
  elsif ($OPT_stages ne "") {
    my @stages = split(",", $OPT_stages);
    my %stages;
    foreach (@stages) {
      die $ERRMSG . "Cannot find subprogram '$_'\n" if (!exists $STAGES{$_});
      $stages{$_} = 1;
    }
    foreach (keys %STAGES) {
      if (exists $stages{$_} && (!defined ${$STAGES{$_}{dependson}} || ${$STAGES{$_}{dependson}} == 1)) {
        $STAGES{$_}{active} = 1;
        $NEED_MAPPINGFILES = 1 if ($STAGES{$_}{needs_mappingfiles});
      }
      else {
        $STAGES{$_}{active} = 0;
      }
    }
  }

}


sub get_mem_in_MB {
  my $mem_str = shift;

  my $mem = $mem_str; chop $mem;
  my $unit = substr($mem_str, length($mem_str)-1);
  if ($unit eq "G") {
    $mem *= 1024;
  }

  return $mem;
}


# not used for now:
sub start_progressbar {
  my ($stage_str, $nr_jobs) = @_;

  my $str = "[MethylScore] " . localtime() . " - $stage_str |";
  $|=1;
  print STDERR "\n" . $str;
  print STDERR "-" x $nr_jobs;
  print STDERR "|\n" . (" " x (length($str)-1)) . "|";
  $|=0;
}


sub print_stage_descr {
  my $stage = shift;
  print STDERR "[MethylScore] " . localtime() . " - STAGE " . $stage . ": " .
    $STAGES{$stage}{"descr"} . "\n";
}


sub pipeline {

  my $s=0;
  for my $stage (sort{$STAGES{$a}{order}<=>$STAGES{$b}{order}} keys %STAGES) {

    next if (!$STAGES{$stage}{active});

    my $first_stage = 0;
    $first_stage = 1 if ($s++ == 0);


    if ($stage eq "dedup") {

        print_stage_descr($stage);

        my @commands;

        foreach my $sample (sort keys %SAMPLES) {
          next if (-e "$MAPPING_FOLDER/$sample/done" && !$CFG::FORCE_RERUN);
          unlink("$MAPPING_FOLDER/$sample/done") if (-e "$MAPPING_FOLDER/$sample/done");

          # check executable first:
          abort_pipeline("Cannot find $CFG::SCRIPT_PATH/merge_and_dedup.sh!\n") if (! -e "$CFG::SCRIPT_PATH/merge_and_dedup.sh");

          #start_progressbar("Processing mapping files", scalar(keys %SAMPLES));

          # pass over 2/3 of stage's maxmem to merge_and_dedup which uses Java:
          my $mem = sprintf("%.0f", get_mem_in_MB($STAGES{"dedup"}{"mem"}) * 2 / 3);
          my $logfile = "$MAPPING_FOLDER/$sample/log";
          mkdir("$MAPPING_FOLDER$sample") if (! -e "$MAPPING_FOLDER$sample");

          my $cmd = "$CFG::SCRIPT_PATH/merge_and_dedup.sh " .
                     $stage . " " .
                     $MAPPING_FOLDER . " " .
                     $sample . " " .
                     join(",", @{$SAMPLES{$sample}{"files"}}) . " " .
                     $SAMPLES{$sample}{"reference"} . " " .
                     $CFG::EXTBIN_PATH . " " .
                     $CFG::FORCE_RERUN . " " . # TODO delete
                     $CFG::REMOVE_INTMED_FILES . " " .
                     $mem . " " .
                     $BAMTOOLS_SPLIT . " " .
                     $CFG::SAMTOOLS . " " .
                     $CFG::STATISTICS . " " .
                     $CFG::SCRIPT_PATH . " " .
                     $CFG::DO_DEDUP . " " .
                     (defined $CFG::ROI ? $CFG::ROI : "");

          push @commands, new_command($cmd, $stage, $logfile);
        }

        run_commands(\@commands);

    }

    elsif ($stage eq "consensus") {

        my %jobs_todo = %SAMPLES;
        while (scalar(keys %jobs_todo) > 0) {

            my @commands;

            foreach my $sample (sort keys %jobs_todo) {

              if (! -e "$MAPPING_FOLDER/$sample/done") {
                if ($first_stage) {
                  abort_pipeline("Stage 'dedup' not finished for sample $sample ($MAPPING_FOLDER$sample/done missing)\n");
                }
                else {
                  sleep 1;
$SLEEP_TOTAL++;
                  check_finished_jobs();
                  next; # wait until previous stage will be finished
                }
              }

              print_stage_descr($stage) if (scalar(keys %jobs_todo) == scalar(keys %SAMPLES));

              my @chr_paths = glob "$MAPPING_FOLDER/$sample/split/*/";
              $SAMPLES{$sample}{"chrs"} = scalar(@chr_paths);
              my $chrs_per_sample = scalar(@chr_paths);
              foreach my $chr_path ( @chr_paths ) {

                chop $chr_path;
                my $chr = substr($chr_path, rindex($chr_path, "/")+1);

                #$UNION_CHROMS{$chr} = () if (! exists $UNION_CHROMS{$chr});
                push @{$UNION_CHROMS{$chr}}, $sample;

                if (-e "$CONSENSUS_FOLDER/$sample/$chr/done" && !$CFG::FORCE_RERUN) {
                    $chrs_per_sample--;
                    next;
                }
                unlink("$CONSENSUS_FOLDER/$sample/$chr/done") if (-e "$CONSENSUS_FOLDER/$sample/$chr/done");

                # check executable first:
                abort_pipeline("Cannot find $CFG::SCRIPT_PATH/consensus.sh!\n") if (! -e "$CFG::SCRIPT_PATH/consensus.sh");

                my $indir  = $MAPPING_FOLDER.$sample."/split/".$chr;
                my $outdir = $CONSENSUS_FOLDER.$sample."/".$chr;
                make_path($outdir) if (! -e $outdir);

                # check if $indir is not empty and make sure indir files are completely copied:
                my $bam_ok=0; my $fa_ok=0;
                for (glob "$indir/*") {
                  $bam_ok=1 if ($_ =~ /.bam$/);
                  $fa_ok=1 if ($_ =~ /.fa$/);
                }
                abort_pipeline("Cannot find .bam or .fa files in consensus input folder $indir. Re-run dedup stage?\n") if (!$bam_ok || !$fa_ok);

                my $logfile = "$outdir/log";
                my $cmd = "$CFG::SCRIPT_PATH/consensus.sh " .
                           $stage . " " .
                           $CFG::EXTBIN_PATH . " " .
                           $sample . " " .
                           $SAMPLES{$sample}{"libtype"} . " " .
                           $CFG::MIN_QUAL . " " .
                           $CFG::IGNORE_LAST_BP . " " .
                           $CFG::IGNORE_FIRST_BP . " " .
                           $indir . " " .
                           $outdir . " " .
                           $CFG::REMOVE_INTMED_FILES . " " .
                           $CFG::SAMTOOLS;

                push @commands, new_command($cmd, $stage, $logfile);
                $chrs_per_sample--;

              }

              delete $jobs_todo{$sample} if ($chrs_per_sample == 0);
            }

            run_commands(\@commands) if (@commands);

            sleep scalar(keys %jobs_todo);
$SLEEP_TOTAL += scalar(keys %jobs_todo);
          }

    }

    elsif ($stage eq "matrix") {

        if ($first_stage || scalar(keys %UNION_CHROMS) == 0) {
          # Fill %UNION_CHROMS and create empty consensus output files if they're missing:
          foreach my $chr (@CHROM_ORDER) {
            my @folders = glob "$CONSENSUS_FOLDER/*/$chr/";
            if (scalar(@folders) > 0) {
              # If we are here, then mapping data for chromosome $chr exists in at least one sample

              # scan all samples and create empty consensus output file if missing
              foreach my $s (sort {$SAMPLES{$a}{"idx"} <=> $SAMPLES{$b}{"idx"}} keys %SAMPLES) {
                if (!-e "$CONSENSUS_FOLDER/$s/$chr/allC.output") {
                  mkdir("$CONSENSUS_FOLDER/$s/$chr") if (!-e "$CONSENSUS_FOLDER/$s/$chr");
                  system("touch $CONSENSUS_FOLDER/$s/$chr/allC.output");
                  warning("[matrix] No mapping data for sample $s and chrom $chr");
                }
                push @{$UNION_CHROMS{$chr}}, $s;
                system("touch $CONSENSUS_FOLDER/$s/$chr/done")
                  if (! -e "$CONSENSUS_FOLDER/$s/$chr/done");
              }

            }
          }
          abort_pipeline("Cannot find consensus files!\n") if (scalar(keys %UNION_CHROMS) == 0);
        }

        my %jobs_todo = %UNION_CHROMS;
        while (scalar(keys %jobs_todo) > 0) {

            my @commands;
            foreach my $chr (sort keys %jobs_todo) {

                my @samples_done = glob "$CONSENSUS_FOLDER/*/$chr/done";
                if (scalar(@samples_done) == scalar(@{$UNION_CHROMS{$chr}})) {

                  print_stage_descr($stage)
                    if (scalar(keys %jobs_todo) == scalar(keys %UNION_CHROMS));

                  if (-e "$MATRIX_FOLDER/$chr/done" && !$CFG::FORCE_RERUN) {
                    delete $jobs_todo{$chr};
                    next;
                  }
                  unlink("$MATRIX_FOLDER/$chr/done") if (-e "$MATRIX_FOLDER/$chr/done");

                  # check executable first:
                  abort_pipeline("Cannot find $CFG::SCRIPT_PATH/generate_matrix.sh!\n")
                    if (! -e "$CFG::SCRIPT_PATH/generate_matrix.sh");

                  ### generate sample file:
                  my $outdir = $MATRIX_FOLDER.$chr;
                  mkdir ($outdir) if (! -e $outdir);

                  my $samplefile = "$outdir/samples.txt";
                  open S, ">$samplefile";
                  foreach my $s (sort{$SAMPLES{$a}{"idx"}<=>$SAMPLES{$b}{"idx"}} keys %SAMPLES) {
                    if (! -e "$CONSENSUS_FOLDER/$s/$chr/allC.output") {
                      # create empty consensus output if it is missing for a sample/chromosome:
                      warning("[matrix] No mapping data for sample $s and chrom $chr");
                      mkdir("$CONSENSUS_FOLDER/$s/$chr") if (! -e "$CONSENSUS_FOLDER/$s/$chr");
                      system("touch $CONSENSUS_FOLDER/$s/$chr/allC.output");
                    }
                    print S $s . "\t" . "$CONSENSUS_FOLDER/$s/$chr/allC.output" . "\n";
                  }
                  close S;
                  ###########################

                  my $logfile = "$outdir/log";
                  my $cmd = "$CFG::SCRIPT_PATH/generate_matrix.sh " .
                             $stage . " " .
                             $CFG::BIN_PATH . " " .
                             $samplefile . " " .
                             $outdir . " " .
                             $chr;

                  push @commands, new_command($cmd, $stage, $logfile);
                  delete $jobs_todo{$chr};
                }
                else {
                    # at least one done file is missing
                    abort_pipeline("Stage 'consensus' not finished for chrom $chr " .
                                   "(done file(s) missing)\n") if ($first_stage);
                    sleep 1;
$SLEEP_TOTAL++;
                    check_finished_jobs();
                }

            } # foreach chr

            run_commands(\@commands) if (@commands);

            sleep scalar(keys %jobs_todo);
$SLEEP_TOTAL += scalar(keys %jobs_todo);
        }
    }

    elsif ($stage eq "matrixWG") {

        wait_until_all_jobs_finished();

        print_stage_descr($stage);

        # SERIALLY: merge chromosomes into global genome matrix:
        if (! -e "$MATRIX_FOLDER/done" || $CFG::FORCE_RERUN) {

          if ($CFG::FORCE_RERUN) {
            unlink("$matrixfile*") if (glob "$matrixfile*" ne "");
            unlink("$MATRIX_FOLDER/done") if (-e "$MATRIX_FOLDER/done");
          }

          # make sure there are all matrices (done files might be present a
          # bit faster than the actual matrix that has to be copied over network/filesystems!)
          my $file_str = "";
          foreach my $chrom (@CHROM_ORDER) {
print DEBUG "chrom $chrom\n" if ($DEBUG);
            next if (! -e "$MATRIX_FOLDER/$chrom");

            my $iter=0;
            while (! -e "$MATRIX_FOLDER/$chrom/genome_matrix.$chrom.tsv" &&
                   ! -e "$MATRIX_FOLDER/$chrom/done") {
print DEBUG "iter $iter\n" if ($DEBUG);
              if (++$iter>10) {
                abort_pipeline("Cannot find $MATRIX_FOLDER/$chrom/genome_matrix.$chrom.tsv!")
              }
              sleep $iter;
            }
            $file_str .= " $MATRIX_FOLDER/$chrom/genome_matrix.$chrom.tsv";
print DEBUG "file_str $file_str\n" if ($DEBUG);
          }

          abort_pipeline("No chromosome specific matrix found!") if ($file_str eq "");
          my $ret = system("cat $file_str | awk '{if (\$0 !~ /^#/ || (\$0 ~ /^#/ && NR==1))
          print}' > $matrixfile 2>> $MATRIX_FOLDER/log");
          abort_pipeline("Cannot create whole-genome matrix! See $MATRIX_FOLDER/log\n")
            if ($ret != 0);

          system("touch $MATRIX_FOLDER/done");
          system("rm $file_str") if ($CFG::REMOVE_INTMED_FILES);

        }

    }

    elsif ($stage eq "MRs") {

        print_stage_descr($stage);

        my @commands;

        foreach my $sample (sort{$SAMPLES{$a}{"idx"}<=>$SAMPLES{$b}{"idx"}} keys %SAMPLES) {
            next if (-e "$MR_FOLDER/$sample/done" && !$CFG::FORCE_RERUN);
            unlink("$MR_FOLDER/$sample/done") if (-e "$MR_FOLDER/$sample/done");

            # check executable first:
            abort_pipeline("Cannot find $CFG::SCRIPT_PATH/call_MRs.sh!\n")
              if (! -e "$CFG::SCRIPT_PATH/call_MRs.sh");

            # check matrix and uncompress if needed:
            if (! -e $matrixfile) {
              if (-e "$matrixfile.gz") {
                my $cmd = $CFG::EXTBIN_PATH . "/bzp -d $matrixfile.gz";
                my $ret = system($cmd);
                if ($ret != 0) {
                  abort_pipeline("bgzip error while decompressing genome matrix\n") if ($ret!=0);
                }
              }
              else {
                abort_pipeline("Cannot find matrix file $matrixfile\n");
              }
            }

            my $logfile = "$MR_FOLDER$sample/log";
            mkdir("$MR_FOLDER$sample") if (! -e "$MR_FOLDER$sample");

            my $cmd = "$CFG::SCRIPT_PATH/call_MRs.sh " .
                       $stage . " " .
                       $CFG::BIN_PATH . " " .
                       $sample . " " .
                       ($SAMPLES{$sample}{"idx"}+1) . " " .
                       $CFG::MIN_COVERAGE . " " .
                       $CFG::DESERT_SIZE . " " .
                       $CFG::MERGE_DIST . " " .
                       (1.0*$CFG::TRIM_METHRATE/100) . " " .
                       $matrixfile . " " .
                       $MR_FOLDER.$sample . " " .
                       $CFG::HUMAN . " " .
                       ($CFG::MR_MIN_C > 0 ? $CFG::MR_MIN_C : -1) . " " .
                       $CFG::MR_PARAMS;

            push @commands, new_command($cmd, $stage, $logfile);
        }

        run_commands(\@commands);

    }

    elsif ($stage eq "igv") {

        print_stage_descr($stage);

        my @commands;

        my @MRfiles = glob "$MR_FOLDER/*/MRs.bed";
        while (scalar(@MRfiles) < scalar(keys %SAMPLES) || (! -e $matrixfile && ! -e "$matrixfile.gz")) {
          if ($first_stage) { abort_pipeline("Input files for stage igv not complete"); }
          check_finished_jobs();
          sleep 1;
$SLEEP_TOTAL++;
          @MRfiles = glob "$MR_FOLDER/*/MRs.bed";
        }

        if (! -e "$IGV_FOLDER/done" || $CFG::FORCE_RERUN) {
          unlink("$IGV_FOLDER/done") if (-e "$IGV_FOLDER/done" && $CFG::FORCE_RERUN);
          mkdir($IGV_FOLDER) if (! -e $IGV_FOLDER);

          # check executable first:
          abort_pipeline("Cannot find $CFG::SCRIPT_PATH/igv.sh!\n") if (! -e "$CFG::SCRIPT_PATH/igv.sh");

          # check matrix and ...
          $matrixfile .= ".gz" if (! -e $matrixfile && -e "$matrixfile.gz");

          # ... uncompress if needed:
          # if ($matrixfile =~ /.gz$/) {
      	  #   if (! -e $matrixfile && -e substr($matrixfile, 0, length($matrixfile)-3)) {
      	  #     $matrixfile = substr($matrixfile, 0, length($matrixfile)-3);
      	  #   }
      	  #   else {
  	      #     abort_pipeline("Cannot find $matrixfile") if (! -e $matrixfile);
          #     my $cmd = $CFG::EXTBIN_PATH . "/bzp -d $matrixfile";
          #     my $ret = system($cmd);
          #     if ($ret != 0) {
          #       abort_pipeline("bgzip error while decompressing genome matrix\n") if ($ret!=0);
          #     }
  	      #     else {
  		    #       $matrixfile = substr($matrixfile, 0, length($matrixfile)-3);
  	      #     }
          #   }
  	      # }

          my $MRfiles = "";
          foreach my $sampleID (sort{$SAMPLES{$a}{"idx"}<=>$SAMPLES{$b}{"idx"}} keys %SAMPLES) {
              abort_pipeline("Cannot find $MR_FOLDER$sampleID/MRs.bed\n") if (! -e "$MR_FOLDER$sampleID/MRs.bed");
              $MRfiles .= "$MR_FOLDER$sampleID/MRs.bed ";
          }


          my $logfile = "$IGV_FOLDER/log";
          my $cmd = "$CFG::SCRIPT_PATH/igv.sh " .
                     $stage . " " .
                     $CFG::SCRIPT_PATH . " " .
                     $matrixfile . " " .
                     $IGV_FOLDER . " " .
                     "\"$MRfiles\" ";

          push @commands, new_command($cmd, $stage, $logfile);
          run_commands(\@commands);
        }

    }

    elsif($stage eq "DMRs") {

        wait_until_all_jobs_finished();

        print_stage_descr($stage);

        if (! -e "$DMR_FOLDER/done" || $CFG::FORCE_RERUN) {
            unlink("$DMR_FOLDER/done") if (-e "$DMR_FOLDER/done" && $CFG::FORCE_RERUN);

            # write sample sheet in a file:
            my $samplefile = "$DMR_FOLDER/samplesheet.tsv";
            open S, ">$samplefile";
            foreach my $sampleID (sort{$SAMPLES{$a}{"idx"}<=>$SAMPLES{$b}{"idx"}} keys %SAMPLES) {

                if ($first_stage && ! -e "$MR_FOLDER$sampleID/MRs.bed") {
                  abort_pipeline("Cannot find MR file of sample $sampleID");
                }

                my $iter=0;
                while (! -e "$MR_FOLDER$sampleID/MRs.bed") {
                  print STDERR "Waiting for $MR_FOLDER$sampleID/MRs.bed\n";
                  sleep ++$iter;
$SLEEP_TOTAL += $iter;
                  abort_pipeline("Waited " . ($iter*$iter/2) . " seconds for MR file of sample $sampleID, ".
                                 "but it's still not there (It could take some time to copy the file through ".
                                 "the network. You can try to run the pipeline or the MR stage again)")
                                 if ($iter > 50);
                }

                print S "$sampleID\t" .
                        ($SAMPLES{$sampleID}{"idx"}+4) . "\t" .
                        "$MR_FOLDER$sampleID/MRs.bed\n";
            }
            close S;


            ##### STEP 1 : split MR file into chunks for cluster parallelization:
            if ( ! -e "$DMR_FOLDER/batches/done" || $CFG::FORCE_RERUN) {
              unlink("$DMR_FOLDER/batches/done") if (-e "$DMR_FOLDER/batches/done" && $CFG::FORCE_RERUN);
              print STDERR "[MethylScore] " . localtime() . " - ... Split MRs into chunks\n" if ($OPT_verbose);

              mkdir("$DMR_FOLDER/batches") if (! -e "$DMR_FOLDER/batches");

              abort_pipeline("Cannot find $CFG::BIN_PATH/split_MRfile\n") if (! -e "$CFG::BIN_PATH/split_MRfile");
              my $cmd = "$CFG::BIN_PATH/split_MRfile " .
                               "$samplefile " .
                               "$DMR_FOLDER/batches/MRbatch " .
                               "$CFG::MR_BATCH_SIZE " .
                               "2>> $DMR_FOLDER/log.split";
              system("echo \"$cmd\n\" >> $DMR_FOLDER/log.split");
              my $ret = system($cmd);

              if ($ret != 0) {
                abort_pipeline("split MR file failed\n");
              }
              else {
                system("touch $DMR_FOLDER/batches/done");
                system("rm $DMR_FOLDER/log.split") if ($CFG::REMOVE_INTMED_FILES);
              }
            }

            my @MR_batches = glob "$DMR_FOLDER/batches/MRbatch*[0-9]";
              abort_pipeline("Cannot find MR batches\n") if (scalar(@MR_batches) == 0);


            ##### STEP 2 : find DMRs IN PARALLEL:
            print STDERR "[MethylScore] " . localtime() . " - ... Calling DMRs\n" if ($OPT_verbose);
            my @commands;
            for my $mr_batch (@MR_batches) {
                foreach my $c (sort keys %OPT_contexts) {
                    if (! -e "$mr_batch.$c.out/done" || $CFG::FORCE_RERUN) {
                        unlink("$mr_batch.$c.out/done")
                          if (-e "$mr_batch.$c.out/done" && $CFG::FORCE_RERUN);

                        # check executable first:
                        abort_pipeline("Cannot find $CFG::BIN_PATH/dmrs-contexts!\n") if (! -e
                          "$CFG::BIN_PATH/dmrs-contexts");
                        # check matrix file:
                        if (! -e $matrixfile) {
                          if (-e "$matrixfile.gz") {
                            $matrixfile .= ".gz";
                          }
                          else {
                            abort_pipeline("Cannot find $matrixfile!\n");
                          }
                        }

                        my $cluster_options = "-i " . $CFG::CLUSTER_MIN_METH_DIFF .
                                             " -j " . $CFG::CLUSTER_MIN_METH . " ";
                        if ($c eq "CG") {
                          $cluster_options = "-i " . $CFG::CLUSTER_MIN_METH_DIFF_CG .
                                            " -j " . $CFG::CLUSTER_MIN_METH_CG . " ";
                        } elsif ($c eq "CHG") {
                          $cluster_options = "-i " . $CFG::CLUSTER_MIN_METH_DIFF_CHG .
                                            " -j " . $CFG::CLUSTER_MIN_METH_CHG . " ";
                        } elsif ($c eq "CHH") {
                          $cluster_options = "-i " . $CFG::CLUSTER_MIN_METH_DIFF_CHH .
                                            " -j " . $CFG::CLUSTER_MIN_METH_CHH . " ";
                        }

                        my $logfile = "$mr_batch.log";
                        my $cmd = "$CFG::BIN_PATH/dmrs-contexts " .
                          "-s $samplefile " .
                          "-r $mr_batch " .
                          "-m $matrixfile " .
                          "-c $c " .
                          "-p " . $CFG::MR_FREQ_CHANGE . " " .
                          $cluster_options .
                          "-v " . $CFG::DMR_MIN_COV . " " .
                          "-n " . $CFG::DMR_MIN_C . " " .
                          "-w " . $CFG::SLIDING_WINDOW_SIZE . " " .
                          "-x " . $CFG::SLIDING_WINDOW_STEP . " " .
                          #"-z " . $CFG::THREADS . " " .
                          "-z 1 " .
                          "-B " . $CFG::BIN_PATH . "/betabin_model " .
                          "-T " . $CFG::EXTBIN_PATH . "/tbx " .
                          "-E " . $CFG::EXTBIN_PATH . "/bzp " .
                          "-K " . $CFG::PYTHON_PATH . " " .
                          "-Y " . $CFG::SCRIPT_PATH . "/pv2qv.py " .
                          "--no-post-process " .
                          "-o $mr_batch.$c.out";
                        $cmd .= " -q" if (!$OPT_verbose);
                        $cmd .= " -g" if ($DEBUG);

                        my $call_cmd = "$CFG::SCRIPT_PATH/call_DMRs.sh $stage \"$cmd\" $mr_batch.$c.out";

                        push @commands, new_command($call_cmd, $stage, $logfile);
                    }
                }
            }

            run_commands(\@commands);

            wait_until_all_jobs_finished();


            ##### STEP 3 : merge DMR files:
            print STDERR "[MethylScore] " . localtime() . " - ... Collect DMRs\n" if ($OPT_verbose);

            foreach my $context (sort keys %OPT_contexts) {
              my $DMRfiles = ();
              for (my $batch_nr = 0; $batch_nr < scalar(@MR_batches); ++$batch_nr) {
                my $batch_file = "$DMR_FOLDER/batches/MRbatch.$batch_nr.$context.out/segments.dif";
                $DMRfiles .= " " . $batch_file if (-e $batch_file);
              }
              if ($DMRfiles eq "") {
                warning("Cannot find any DMR files (context $context)!\n");
                system("touch $DMR_FOLDER/batch_empty");
                $DMRfiles .= "$DMR_FOLDER/batch_empty";
              }

              my $cmd = "cat $DMRfiles > $DMR_FOLDER/segments.$context.dif";
              my $ret = system($cmd);
              abort_pipeline("DMR calling (catting segments.$context.dif's) failed\n") if ($ret);

              abort_pipeline("Cannot find $CFG::BIN_PATH/merge_DMRs-contexts\n") if (!-e
                  "$CFG::BIN_PATH/merge_DMRs-contexts");
              my $cmd = "$CFG::BIN_PATH/merge_DMRs-contexts " .
                  "$samplefile " .
                  "$DMR_FOLDER/segments.$context.dif " .
                  "$DMR_FOLDER " .
                  "$OPT_config_file " .
                  "$context " .
                  "2>> $DMR_FOLDER/log.merge";
              system("echo \"$cmd\n\" >> $DMR_FOLDER/log.merge");
              my $ret = system($cmd);

              if ($ret) {
                abort_pipeline("DMR calling (collect step) failed\n");
              }
            }

            # on success:
            remove_tree("$DMR_FOLDER/batches") if ($CFG::REMOVE_INTMED_FILES);
            system("rm $DMR_FOLDER/log.merge") if ($CFG::REMOVE_INTMED_FILES);
            unlink("$DMR_FOLDER/batch_empty") if (-e "$DMR_FOLDER/batch_empty");
            system("touch $DMR_FOLDER/done");

         } # if (run this stage)

     } # end of stage DMR


     # start next stage only if threads are available (i.e. if there are less jobs running or pending than $CFG::THREADS):
     while (scalar(keys %JOBS_PENDING) >= $CFG::THREADS) {

       check_finished_jobs(); # reduces %JOBS_PENDING entries

       sleep 1;
       sleep scalar(keys %JOBS_PENDING) - $CFG::THREADS if (scalar(keys %JOBS_PENDING) >=
         $CFG::THREADS);
$SLEEP_TOTAL += (scalar(keys %JOBS_PENDING) >= $CFG::THREADS) ? scalar(keys %JOBS_PENDING) -
         $CFG::THREADS + 1 : 1;
     }

  } # for each stage


  wait_until_all_jobs_finished();

  $Thread_pool->shutdown if (!$CFG::CLUSTER);
  if (-e $TMP_FOLDER) { remove_tree("$TMP_FOLDER"); }

  print STDERR "\n[MethylScore] " . localtime() . " - FINISHED SUCCESSFULLY\n" if ($OPT_verbose);

if ($DEBUG) { print DEBUG "\nTotal sleeping time: $SLEEP_TOTAL\n"; }

}




sub wait_until_all_jobs_finished {
  my $sec = 0;
  while (scalar(keys %JOBS_PENDING) > 0) {
    check_finished_jobs();
    sleep scalar(keys %JOBS_PENDING);
$SLEEP_TOTAL += scalar(keys %JOBS_PENDING);
    $sec += scalar(keys %JOBS_PENDING);
  }
print DEBUG "waited $sec secs for all jobs to finish\n" if ($DEBUG);
}



sub new_command {
  my ($str, $stage, $logfile, $is_binary) = @_;

  my %cmd;
  $cmd{"str"} = $str;
  #$cmd{"cpus"} = $STAGES{$stage}{"cpus"};
  $cmd{"time"} = $STAGES{$stage}{"time"};
  $cmd{"mem"} = $STAGES{$stage}{"mem"};
  $cmd{"progressbar"} = $STAGES{$stage}{"progressbar"};
  $cmd{"id"} = $JOBID++;
  $cmd{"binary"} = $is_binary ? 1 : 0;
  $cmd{"stage"} = $stage;
  $cmd{"logfile"} = $logfile;

  return \%cmd;
}



sub run_commands {
  my ($c) = @_;
  my @commands = @{$c};

  if ($CFG::CLUSTER) {

print DEBUG "[run] batch of " . @commands . " coming in\n" if ($DEBUG);

      while (scalar(@commands) > 0) {
        my $running_jobs = parse_submitted_jobs();
        my $nr_jobs_to_submit = min($CFG::THREADS-$running_jobs, scalar(@commands));

        # submit ($CFG::THREADS - $running_jobs) new jobs:
        for (my $i=0; $i<$nr_jobs_to_submit; ++$i) {

          my %command = %{$commands[0]};

          # remove logfile if CFG::FORCE_RERUN
          unlink $command{"logfile"} if (-e $command{"logfile"});

          my $cmd_to_submit;
          if (!$SGE) {
            # run on GMI cluster:
            $cmd_to_submit = "qsub " .
                             "-l mem=" . $command{"mem"} . " " .
                             "-l walltime=" . $command{"time"} . " " .
                             "-P " . $CFG::CLUSTER_PROJECT . " " .
                              "-o " . $command{"logfile"} . " " .
                              "-j oe " .
                              "-V " .
                              "-- " .
                              #($command{"binary"} ? " -- " : "") .
                              $command{"str"};
          }
          else {
            # run on default SGE cluster:
            $cmd_to_submit = "qsub " .
                             "-o " . $command{"logfile"} . " " .
                             "-j y " .
                             "-l h_vmem=" . $command{"mem"} . " " .
                             #"-pe smp " . $command{"cpus"} . " " .
                             ($command{"binary"} ? "-b y " : "-b n ") .
                             $command{"str"};
          }

if ($DEBUG) { print DEBUG "$cmd_to_submit\n"; }
          # write command into logfile:
          system("echo \"$cmd_to_submit\n\" >> " . $command{"logfile"});

          # submit command and fetch stdout:
          my $outstr = `$cmd_to_submit`;
if ($DEBUG) { chop $outstr; print DEBUG "[run] SUBMITTED:   $cmd_to_submit ($outstr with ".$command{"mem"}."MB RAM)\n"; }

          # parse jobid from stdout:
          my $jobid;
          if (!$SGE) {
            chomp($jobid = $outstr);
          }
          else {
            $outstr =~ /Your job (\d+) \(/;
            $jobid = $1;
          }

          # write jobid into logfile:
          system("echo \"\nJOBID: $jobid\" >> " . $command{"logfile"});

print DEBUG "[run] jobid: $jobid\n" if ($DEBUG);

          # log which command was submitted:
          $JOBS_PENDING{$jobid} = \%command;
print DEBUG "[run] jobs pending: ".scalar(keys %JOBS_PENDING)."\n" if ($DEBUG);

          # delete submitted command from @commands that still have to be submitted:
          shift @commands;
        }

        sleep 1+$running_jobs;
$SLEEP_TOTAL += 1+$running_jobs;
      }

  }
  else {
      ### run locally:
      # max no. threads is taken care of by Thread::Pool
      foreach my $cmd (@commands) {
        my %command = %{$cmd};
        system("echo -e \"" . $command{"str"} . "\n\" >> " . $command{"logfile"});
        my $str = $command{"str"} . " 2>> " . $command{"logfile"} . " >> " . $command{"logfile"};
        my $jobid = $Thread_pool->job("$str");
        $JOBS_PENDING{$jobid} = \%command;
print DEBUG "[run] RUNNING Job $jobid: $str\n" if ($DEBUG);
print DEBUG "[run] jobs pending: ".scalar(keys %JOBS_PENDING)."\n" if ($DEBUG);
      }
  }

}



###########################################################
### Parse qstat, get number of submitted jobs (running or
### the queue)
sub parse_submitted_jobs {
print DEBUG "[parse] jobs pending: ".scalar(keys %JOBS_PENDING)."\n" if ($DEBUG);
  %JOBS_RUNNING = ();

  my $qstat_output = `qstat -u \$USER`;
print DEBUG "$qstat_output--------------\n" if ($DEBUG);
  my @lines = split"\n", $qstat_output;

  my $nr_jobs = 0;

  for my $line (@lines) {
    my @fields = split" ", $line;
    next if ($line =~ /^[jJ-]/); # qstat starts with job-ID (1st line) and ---(2nd line)

    my $jobid = $fields[0];
    $jobid .= $CFG::JOBID_TAIL if (!$SGE);
#print DEBUG "$line (JOBID: $jobid, jobs_pending: ".scalar(keys %JOBS_PENDING).")\n";
    my $status = $fields[4];

    if (exists $JOBS_PENDING{$jobid}) {
      if ($SGE && $status =~ /E/) {
        # Job failed:
        abort_pipeline("Job $jobid failed:\n" . $JOBS_PENDING{$jobid}{"str"});
      }
      $nr_jobs++;
      $JOBS_RUNNING{$jobid} = undef; # just to create the hash entry
    }
  }

print DEBUG "[parse] jobs running: $nr_jobs ".scalar(keys %JOBS_RUNNING)."\n" if ($DEBUG);
  return $nr_jobs;
}



sub check_finished_jobs {
print DEBUG "check jobs, pending: " . scalar(keys %JOBS_PENDING) . "\n" if ($DEBUG);

  if ($CFG::CLUSTER) {
    parse_submitted_jobs(); # sets %JOBS_RUNNING
  }

  foreach my $jobid (keys %JOBS_PENDING) {
    if ($CFG::CLUSTER) {
      next if (exists $JOBS_RUNNING{$jobid});

      if ($SGE) {
        # update JOBS_FINISHED if jobs missing. If not, do not parse
        parse_accounting_file() if (!exists $JOBS_FINISHED{$jobid});
      }
      else {
        parse_qstat_xf($jobid);
      }
    }

    my $job_status = get_job_status($jobid);
    # undef: job still running
    # 0    : Job successful
    # !=0  : Job faield

print DEBUG "  [check] job $jobid status: $job_status\n" if ($CFG::CLUSTER && $DEBUG);

  #print "Job $jid has exit status $job_status\n";

    next if (! defined $job_status);

    # write max memory usage and runtime of job into job log file:
    if ($CFG::CLUSTER) {
      open F, ">>", $JOBS_PENDING{$jobid}{"logfile"};
      print F "\n\nMem requested: " . $JOBS_PENDING{$jobid}{"mem"} . "\n";
      #print F "CPUs requested: " . $JOBS_PENDING{$jobid}{"cpus"} . "\n";
      print F "Max mem (Bytes): " . $JOBS_PENDING{$jobid}{"maxmem"} . "\n";
      print F "Runtime (sec): " . $JOBS_PENDING{$jobid}{"runtime"} . "\n";
      #print F "Logfile: " . $JOBS_PENDING{$jobid}{"logfile"} . "\n";
      close F;
    }

    if ($job_status != 0) { # Job failed
#      if (!$CFG::CLUSTER) {
#        # when run locally, the corresponding script will output an error message
#        abort_pipeline("Command failed:\n".$JOBS_PENDING{$jobid}{"str"});
#      }
#      else {
        # on cluster, output last 5 lines of the job's logfile:
        my $logfile = $JOBS_PENDING{$jobid}{"logfile"};
        abort_pipeline("Command failed:\n".$JOBS_PENDING{$jobid}{"str"}."\n\n" .
                       "Check out the logfile: $logfile\n" .
                       "Tail of logfile:\n'''''\n" . `tail $logfile` .
                       "'''''");
#      }
    }
    else { # Job finished successfully, delete it from @JOBS_PENDING
      delete $JOBS_PENDING{$jobid};
      delete $JOBS_FINISHED{$jobid};
      delete $JOBS_RUNNING{$jobid} if (exists $JOBS_RUNNING{$jobid});
    }
  }
print DEBUG "  [check] ... now still pending: ".scalar(keys %JOBS_PENDING)."\n" if ($DEBUG);
}


# parses 'qstat -x -f $jobid' on GMI cluster to get the exit status and runtime/mem stats:
sub parse_qstat_xf {
  my ($jobid) = @_;

  my $qstat_output = `qstat -x -f $jobid`;
if ($DEBUG) { print DEBUG "qstat -x -f $jobid:\n------------------------\n" .
"$qstat_output\n------------------------\n"; }
  abort_pipeline("Cannot retrieve qstat information: qstat -x -f $jobid") if (!$qstat_output);

  if ($qstat_output =~ /Exit_status = (\d+)/) {
    $JOBS_FINISHED{$jobid} = $1;
  } else {
    $JOBS_FINISHED{$jobid} = undef;
    return;
  }

  $JOBS_PENDING{$jobid}{"maxmem"} = $1 if ($qstat_output =~ /resources_used.mem = (.*$)/);
  $JOBS_PENDING{$jobid}{"runtime"} = $1 if ($qstat_output =~ /resources_used.walltime = (.*$)/);

print DEBUG "[qstat_xf] found exit code of Job $jobid: ".$JOBS_FINISHED{$jobid}."\n"
  if ($1 && $DEBUG);
}


# return value:
#   undef: job still running
#       0: job successful
#    != 0: job failed
sub get_job_status {
  my ($job_id) = @_;

  if (!$CFG::CLUSTER) {
    return $Thread_pool->result_dontwait($job_id);
  }
  else {
    return $JOBS_FINISHED{$job_id};
  }
}


### parses cluster accounting file for jobs that
### are in JOBS_PENDING but not in JOBS_RUNNING and not in JOBS_FINISHED
sub parse_accounting_file {

    my $jobs_to_find = scalar(keys %JOBS_PENDING) - scalar(keys %JOBS_RUNNING) - scalar(keys %JOBS_FINISHED);
print DEBUG "[acc] Start parsing acc file and look for $jobs_to_find Jobs\n" if ($DEBUG);
    abort_pipeline("No. jobs to find is 0 in parse_accounting_file") if ($jobs_to_find == 0);

    my $tail_iterations = 0;
    my $tail_len = $ACCOUNTING_FILE_TAIL_LEN;
    do {
      my $accfile = `tail -n $tail_len $ACCOUNTING_FILE`;

      my @lines = split"\n", $accfile;
print DEBUG "[acc] no lines of acc file: ".scalar(@lines)."\n" if ($DEBUG);

      my $f=0;
      foreach my $line (@lines) {

        next if ($line =~ /^[#\s]/);

        my @fields = split":", $line;
        my $jobid = $fields[5];
print DEBUG "[acc] 5th line: ".$line." jobid: $jobid\n" if ($f++ == 5 && $DEBUG);

        if (exists $JOBS_RUNNING{$jobid} || exists $JOBS_FINISHED{$jobid}) {
          $JOBS_FINISHED{$jobid} = undef if (exists $JOBS_RUNNING{$jobid});
          next;
        }

        if (exists $JOBS_PENDING{$jobid}) {
          my $exit_code = $fields[12];
          $JOBS_PENDING{$jobid}{"maxmem"} = int(sprintf("%.0f", $fields[42])); # bytes
          $JOBS_PENDING{$jobid}{"runtime"} = int(sprintf("%.0f", $fields[14]));

          $JOBS_FINISHED{$jobid} = $exit_code;
print DEBUG "[acc] found Job $jobid in accounting file, exit code: $exit_code\n" if ($DEBUG);

          $jobs_to_find--;
print DEBUG "[acc] finished parsing accounting file\n" if ($jobs_to_find==0 && $DEBUG);
          return if ($jobs_to_find == 0);
        }

      }

      sleep $ACCOUNTING_FILE_SLEEP;
$SLEEP_TOTAL+=$ACCOUNTING_FILE_SLEEP;
      $tail_len *= 2;
      $tail_iterations++;
    } while ($tail_len < $MAX_ACCFILE_TAIL_LEN);

    abort_pipeline("Tried $tail_iterations times to parse cluster accounting file, " .
                   "but could not find all submitted jobs\n");
}


sub warning {
  my ($war) = @_;

  print WARNINGS $WARNMSG . $war . "\n";
  $WARNINGS = 1;
}



sub abort_pipeline {
  my ($err) = @_;
  if ($CFG::CLUSTER) {
    # kill all currently submitted jobs:
    parse_submitted_jobs();
    foreach my $jobid (keys %JOBS_RUNNING) {
      system("qdel $jobid");
    }
  }
  else {
    $Thread_pool->abort if (defined $Thread_pool);
  }

  if (-e $TMP_FOLDER) { remove_tree("$TMP_FOLDER"); }

  print STDERR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" .
               "ERROR!\n$err\n" .
               "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" .
               "MethylScore finished with ERRORS\n";
  print STDERR "...and there are warnings ($CFG::PROJECT_FOLDER/warnings.txt)\n" if ($WARNINGS);
  exit 1;
}


### Read command line parameters --------------------------------------------------------------
sub GetCom {

  my $subprogs = ""; my $i=0;
  for my $subprog (sort{$STAGES{$a}{order}<=>$STAGES{$b}{order}} keys %STAGES) {
    $subprogs .= "\n                                " if ($i>0 && $i % 4 == 0);
    $subprogs .= $subprog . " ";
    $i++;
  }

	my @usage = ("$0

MethylScore - Detection of differential methylation using NGS
version $version

Mandatory parameter:
-c|--config         <path>   Path to config file

Other parameters (NOTE: they override config file parameters):
-a|--stages         STRING   (comma-separated) Subprograms to run, choose from:
                                $subprogs
                             'all' selects all subprograms
                             default: all but 'igv'
-p|--project        <path>   Project folder (will be created if not existant)
-s|--sample-sheet   <path>   Path to sample sheet file

-t|--threads        INT      Maximum threads/CPUs to use
-l|--local                   Run locally and not on cluster system
-f|--force-rerun             Force to re-run already performed steps
-d|--donot-dedup             Skip read de-duplication
-k|--keep-files              Do not delete intermediate files
-q|--quiet                   Quiet mode
-h|--help                    Prints this help message

See manual for further explanations.

Copyright (C) 2016-2020, Computomics GmbH
\n");

  if (@ARGV == 0) {
    print STDERR (@usage);
    exit 1;
  }

	my $command = $0 . " " . join(" ", @ARGV);

  GetOptions(
    "config=s" => \$OPT_config_file,
    "s|sample-sheet=s" => \$OPT_sample_sheet,
    "project=s" => \$OPT_project_folder,
    "a|stages=s" => \$OPT_stages,
    "threads=i" => \$OPT_threads,
    "local" => \$OPT_run_local,
    "force-rerun" => \$OPT_force_rerun,
    "d|donot-dedup" => \$OPT_donot_dedup,
    "keep-files" => sub { $OPT_rm_intmd = 0; },
    "quiet" => sub{ $OPT_verbose = 0; },
    "sge" => \$SGE,
    "debug" => \$DEBUG,
    "help" => sub { print STDERR (@usage); exit 0; }
  );

  abort_pipeline("Please specify config file!") unless defined $OPT_config_file;

#  if ($OPT_project_folder) {
#    $OPT_project_folder = abs_path($OPT_project_folder) or
#      abort_pipeline("Cannot create project folder");
#  }

  return $command;
}
