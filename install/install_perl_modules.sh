#!/bin/bash
set -e

OPTIND=1 # Reset getopts
SHELL_FILE=~/.bash_profile
PERL_BIN=perl

while getopts "h?s:p:" opt; do
	case "$opt" in
	h|\?)
		printf "\nconfigure options:\n-s   shell config file (default: ~/.bash_profile)\n-p   perl executable (default: perl)\n\n"
		exit 0
		;;
	s)	SHELL_FILE=$OPTARG
		;;
	p)	PERL_BIN=$OPTARG
		;;
	esac
done

shift $((OPTIND-1))
[ "$1" = "--" ] && shift


perl -v

# where to install perl modules:
export LOCAL_PERL_PATH=~/MethylScore/lib/perl5
if [ ! -e $LOCAL_PERL_PATH ]; then mkdir -p $LOCAL_PERL_PATH; fi

wget -O - http://cpanmin.us | perl - -l $LOCAL_PERL_PATH App::cpanminus local::lib
eval `perl -I $LOCAL_PERL_PATH/lib/perl5 -Mlocal::lib=$LOCAL_PERL_PATH`

if ! grep -q "eval \`perl -I $LOCAL_PERL_PATH/lib/perl5" $SHELL_FILE; then
  echo "eval \`perl -I $LOCAL_PERL_PATH/lib/perl5 -Mlocal::lib=$LOCAL_PERL_PATH\`" >> $SHELL_FILE
  echo "export MANPATH=$LOCAL_PERL_PATH/man:\$MANPATH" >> $SHELL_FILE
  source $SHELL_FILE
fi


# all modules that have to be installed:
modules="Getopt::Long Thread::Pool Thread::Conveyor::Array Thread::Conveyor::Throttled Thread::Conveyor::Tied Thread::Tie::Array Config::Simple File::Path File::Copy File::Tee File::Which Algorithm::Cluster Data::Dumper List::Util Compress::BGZF::Reader PAR::Packer experimental"

$PERL_BIN -e 'use ExtUtils::Installed; my $inst = ExtUtils::Installed->new(); my @modules = $inst->modules(); print join(" ", @modules) . "\n";' > tmp_installedPerlModls


# install missing modules:
for m in $modules; do
  if ! grep $m tmp_installedPerlModls; then
    cpanm --force $m
  fi
done

source $SHELL_FILE

rm tmp_installedPerlModls
