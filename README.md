![Integration test](https://github.com/Computomics/MethylScore/workflows/Integration%20test/badge.svg?branch=nextflow)
![Docker Container](https://github.com/Computomics/MethylScore/workflows/Build%20Container/badge.svg?branch=nextflow)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)
# MethylScore-nf
A [nextflow](https://www.nextflow.io/) implementation of Computomics' MethylScore pipeline

# Usage

## Setting up the environment

The only software dependencies that are required to run MethylScore are [nextflow](https://www.nextflow.io/), and one of the following: [Singularity](https://www.sylabs.io/singularity/), [Docker](https://www.docker.com/), [Charliecloud](https://hpc.github.io/charliecloud/) or [Podman](https://podman.io/).

## Running the pipeline

To run the pipeline, you have to provide atleast `--SAMPLE_SHEET` and `--GENOME` (same as the one which you mapped the reads against)

```bash
nextflow run Computomics/MethylScore --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa -profile cbe,singularity
```

## Parameters

Pipeline parameters are initialized to their [default values](https://github.com/Computomics/MethylScore/raw/nextflow/example_config.yaml).

Individual parameters can be passed to the pipeline on the commandline. For example if you want methylated regions to be visualized as [IGV](https://software.broadinstitute.org/software/igv/) tracks:

```bash
nextflow run Computomics/MethylScore --SAMPLE_SHEET=/path/to/samplesheet.tsv --IGV -profile cbe,singularity
```

Alternatively, the repository also contains an [example_config.yaml](https://github.com/Computomics/MethylScore/raw/nextflow/example_config.yaml), which can be used to pass custom parameters to the pipeline using the `-params-file` flag.

```bash
nextflow run Computomics/MethylScore -params-file=/path/to/config.yaml
```

## Profiles
There are three profiles which use different methods to (automatically) resolve external pipeline dependencies.
Note that you have to provide one of them with the `-profile` flag.

### Singularity `-profile singularity` or Charliecloud `-profile charliecloud`
The Singularity and Charliecloud profiles uses a [prebuilt container](https://quay.io/repository/beckerlab/methylscore) that includes all dependencies needed to run the pipeline.

### Docker or Podman `-profile docker` or `-profile podman`
While these are very similar to Singularity or Charliecloud, they are usually not available in Cluster environments.
They are nevertheless very useful for local runs of the pipeline.

Note: with `-profile local,docker` `-profile local,podman` you can easily run the pipeline on your local computer, given you have nextflow and Docker or Podman installed.

### Cluster specific profiles

#### CBE `-profile cbe`
This profile is optimized for the CLIP Batch Environment (CBE) at the Vienna BioCenter and uses `singularity` by default. Enable with `-profile cbe`.

