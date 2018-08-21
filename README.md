# MethylScore-nf
A nextflow implementation of Computomics' MethylScore pipeline

# Usage

```bash
git clone https://github.com/Gregor-Mendel-Institute/MethylScore-nf

module load Nextflow/0.31.1-Java-10.0.2

nextflow /path/to/MethylScore.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --GENOME=/path/to/reference_genome.fa
```
Currently, the reference genome defaults to TAIR9. Please change accordingly if working with other organisms.

Pipeline parameters are included in the nextflow script itself and are initialized to their default values.
Currently, passing a config file is not supported, but individual parameters can be passed to the pipeline.

For example if you want methylated regions to be visualized as igv:

```bash
nextflow /path/to/MethylScore.nf --SAMPLE_SHEET=/path/to/samplesheet.tsv --IGV=true
```
