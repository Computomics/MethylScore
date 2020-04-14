[![Build Status](https://travis-ci.com/Gregor-Mendel-Institute/MethylScore-nf.svg?token=RozNRzpisdeKpeAjRY7S&branch=master)](https://travis-ci.com/Gregor-Mendel-Institute/MethylScore-nf)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Container available](https://img.shields.io/docker/automated/beckerlab/methylscore.svg)](https://hub.docker.com/r/beckerlab/methylscore/)

# MethylScore-nf
A [nextflow](https://www.nextflow.io/) implementation of Computomics' MethylScore pipeline

# Usage

## Setting up the environment

The only software dependencies that are required to run MethylScore-nf are [nextflow](https://www.nextflow.io/), and one of the following: [Singularity](https://www.sylabs.io/singularity/), [Conda](https://conda.io) or [Docker](https://www.docker.com/).

For example if you want to run the pipeline on the CBE cluster using Singularity, do the following on the login node:

```bash
module load nextflow/20.01.0 && module load singularity/3.4.1
```

## Running the pipeline

To run the pipeline, you have to provide atleast `--SAMPLE_SHEET` and `--GENOME` (same as the one which you mapped the reads against)

As long as the pipeline repository is private, you have to use the `-user` flag to authenticate with your Github credentials.

```bash
nextflow run Gregor-Mendel-Institute/MethylScore-nf -user your_github_username --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa -profile cbe,singularity
```

## Parameters

Pipeline parameters are initialized to their [default values](https://github.com/Gregor-Mendel-Institute/MethylScore-nf/raw/master/example_config.yaml).

Individual parameters can be passed to the pipeline on the commandline. For example if you want methylated regions to be visualized as [IGV](https://software.broadinstitute.org/software/igv/) tracks:

```bash
nextflow run Gregor-Mendel-Institute/MethylScore-nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --IGV -profile cbe,singularity
```

Alternatively, the repository also contains an [example_config.yaml](https://github.com/Gregor-Mendel-Institute/MethylScore-nf/raw/master/example_config.yaml), which can be used to pass custom parameters to the pipeline using the `-params-file` flag.

```bash
nextflow run Gregor-Mendel-Institute/MethylScore-nf -params-file=/path/to/config.yaml
```

## Profiles
There are three profiles which use different methods to (automatically) resolve external pipeline dependencies.
Note that you have to provide one of them with the `-profile` flag, in addition to the cluster specific profile.

### Singularity
The Singularity profile uses a [prebuilt docker container](https://hub.docker.com/r/beckerlab/methylscore/) that includes all dependencies needed to run the pipeline.
To run on CBE, you have to `module load singularity/3.4.1` first and then specify `-profile cbe,singularity`

### Conda
The Conda profile uses [Bioconda](https://bioconda.github.io/) to resolve dependencies. It creates a virtual environment and resolves dependencies according to the `environment.yaml` file.
To run on CBE, you have to `module load anaconda3/2019.10` first and then specify `-profile cbe,conda`

### Docker
While this is very similar to the Singularity profile, it is usually not available in Cluster environments.
It is nevertheless very useful for local runs of the pipeline.

Note: with `-profile local,conda` or `-profile local,docker` you can easily run the pipeline on your local computer, given you have nextflow and Conda or Docker installed.
