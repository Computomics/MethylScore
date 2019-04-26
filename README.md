[![Build Status](https://travis-ci.com/Gregor-Mendel-Institute/MethylScore-nf.svg?token=RozNRzpisdeKpeAjRY7S&branch=master)](https://travis-ci.com/Gregor-Mendel-Institute/MethylScore-nf)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.01.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Container available](https://img.shields.io/docker/automated/beckerlab/methylscore.svg)](https://hub.docker.com/r/beckerlab/methylscore/)

# MethylScore-nf
A nextflow implementation of Computomics' MethylScore pipeline

# Usage

## Setting up the pipeline

```bash
git clone https://github.com/Gregor-Mendel-Institute/MethylScore-nf

module load Nextflow
```

Instead of cloning the repository you can let nextflow take care of it:
```bash
nextflow pull Gregor-Mendel-Institute/MethylScore-nf -user your.github.username
```
As the pipeline repository is private, you have to be in the @BeckerLab github team and use the `-user` flag to authenticate.


## Running the pipeline
To run the pipeline on MENDEL you have to provide atleast `--SAMPLE_SHEET` and `--GENOME` (same as the one which you mapped the reads against)

```bash
nextflow MethylScore/main.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa -profile mendel,singularity
```
or if running directly from github:
```bash
nextflow run Gregor-Mendel-Institute/MethylScore-nf -user your.github.username --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa -profile mendel,singularity
```

## Parameters
Pipeline parameters are included in the nextflow script itself and are initialized to their default values.

Individual parameters can be passed to the pipeline on the commandline. For example if you want methylated regions to be visualized as igv tracks:

```bash
nextflow MethylScore/main.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --IGV -profile mendel,singularity
```

The repository also contains `example_config.yaml` which can be edited and then passed to the pipeline using the `-params-file` flag

```bash
nextflow MethylScore/main.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv -params-file=/path/to/config.yaml
```

## Profiles
There are two profiles which use different methods to resolve external pipeline dependencies.
Note that you have to provide one of them with the `-profile` flag.


### singularity
The singularity profile uses a container with a basic linux system that includes all dependencies needed to run the pipeline.

_Mendel only:_
A prebuilt image is located at `/lustre/scratch/projects/becker_common/singularity_images/MethylScore.simg`and is used automatically when `-profile mendel,singularity` is set.


### conda
The conda profile uses the conda system to resolve dependencies. It creates virtual environment and resolves dependencies according to the `environment.yaml` file.
To run on mendel, you have to `module load Miniconda3` first.

Note: with `-profile local,conda` or `-profile local,docker` you can easily run the pipeline on your local computer, given you have nextflow and conda or docker installed.
