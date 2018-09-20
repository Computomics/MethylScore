# MethylScore-nf
A nextflow implementation of Computomics' MethylScore pipeline

# Usage

## Setting up the pipeline

```bash
git clone https://github.com/Gregor-Mendel-Institute/MethylScore-nf

module load Nextflow/0.31.1-Java-10.0.2
```

Instead of cloning the repository you can let nextflow take care of it:
```bash
nextflow pull Gregor-Mendel-Institute/MethylScore-nf -user your.github.username
```
As the pipeline repository is private, you have to be in the @BeckerLab github team and use the `-user` flag to authenticate.


## Running the pipeline
To run the pipeline you have to provide atleast `--SAMPLE_SHEET` and `--GENOME` (same as the one which you mapped the reads against)

```bash
nextflow /path/to/main.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa
```
or if running directly from github:
```bash
nextflow run Gregor-Mendel-Institute/MethylScore-nf -user your.github.username --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa
```

## Parameters
Pipeline parameters are included in the nextflow script itself and are initialized to their default values.

Individual parameters can be passed to the pipeline on the commandline. For example if you want methylated regions to be visualized as igv:

```bash
nextflow /path/to/MethylScore.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --IGV=true
```

The repository also contains `example_config.yaml` which can be edited and then passed to the pipeline using the `-params-file` flag

```bash
nextflow /path/to/MethylScore.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv -params-file=/path/to/config.yaml
```

## Profiles
There are two profiles which use different methods to resolve external pipeline dependencies.
By default, the pipeline uses the singularity profile described below.

### singularity
The singularity profile uses a container with a basic linux system that includes all dependencies needed to run the pipeline.
Because the `Singularity` deployment on MENDEL lacks some features (caused by the compute nodes running a very outdated linux kernel), the image has to be build externally using the `Singularity` and `environment.yaml` files.
```bash
singularity build MethylScore.simg Singularity
```
A prebuilt image is located at `/lustre/scratch/projects/becker_common/singularity_images/MethylScore.simg`and is used automatically when `-profile singularity` is set.

### conda
The conda profile uses the conda system to resolve dependencies. It creates virtual environment and resolves dependencies according to the `environment.yaml` file.
To run on mendel, you have to `module load Miniconda3` first.

Note: with `-profile conda` you can easily run the pipeline on your local computer, given you have nextflow and conda installed and set `executor = 'local'` in the config file.
