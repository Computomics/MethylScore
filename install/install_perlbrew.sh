#!/bin/bash
set -e

usage="Usage: $0\n  <install dir for perlbrew, ABSOLUTE path (e.g. ~/perl5/perlbrew)>\n  <bash shell configuration file (e.g. ~/.bash_profile)>"

if [ -z "$1" -o -z "$2" ]; then echo -e $usage; exit 1; fi
if [[ $1 != /* ]] || [[ $2 != /* ]]; then echo -e "Please use absolute paths!"; exit 1; fi

# install perlbrew to update perl to >=5.14 (for some modules, e.g. Thread::Pool)
export PERLBREW_ROOT=$1         # install directory
BASH_FILE=$2
wget -O - https://install.perlbrew.pl | bash
# alternative:
#curl -L https://install.perlbrew.pl | bash

$PERLBREW_ROOT/bin/perlbrew init
source $PERLBREW_ROOT/etc/bashrc
echo -e "\nsource $PERLBREW_ROOT/etc/bashrc" >> $BASH_FILE

# installation of perl versions requires gcc:
# sudo yum install gcc.x86_64

# continue perl update
$PERLBREW_ROOT/bin/perlbrew install --notest --thread --multi perl-5.22.3
$PERLBREW_ROOT/bin/perlbrew switch perl-5.22.3
source $PERLBREW_ROOT/etc/bashrc

