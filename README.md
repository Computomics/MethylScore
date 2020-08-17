# MethylScore

version 0.1.17

---

Pipeline to call differential methylation simultaneously between multiple samples using next-gen sequencing data.

## Installation

Requirements:
- perl >= 5.14 (for older versions, see section [Older perl versions](#Older-perl-versions))
- perl packages: `bash install/install_perl_packages.sh`
- python >= 3
  - packages scipy and argparse

### Binaries
There are pre-compiled binaries for Ubuntu bionic 18.04 in the folder `bin/`.

### From Source
For now, perl scripts have to be compiled (`src/*pl`), and a C++ program for the MR detection has to be compiled (`src/hmm_mrs`).

You need to have installed:
- perl package `PAR::Packer` (it is included in the install script `install/install_perl_packages.sh`)
- g++

Compile perl scripts and C++ program by:
- executing `compile.sh` in the main folder

### Older perl versions

For perl versions < 5.14:
- install perlbrew by executing `install/install_perlbrew.sh` 
- set the flag `use_perlbrew` in the MethylScore binary file `MethylScore` to 1, so that line 4 looks like:<br>line 4: `use_perlbrew=1`<br> (for `MethylScore-contexts` accordingly)

## Running

`./MethylScore` in the main folder contains an overview of the parameters. The manual (`docs/MethylScore_manual.pdf`) specifies more details and all parameters that can be listed in a config file.

## Manual

See in `docs/MethylScore_manual.pdf`

TODO: Manual does not yet include sequence context specific calling of DMRs.


## LICENSE

Licensed under the Non-Profit Open Software License version 3.0.

If you are a for-profit organization, please contact us to obtain a license at: info@computomics.com

See the file `LICENSE` for the complete license.