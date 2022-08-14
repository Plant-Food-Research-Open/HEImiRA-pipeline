# HEImiRA Pipeline
_Host environment intersection miRNA annotation_ pipeline

## Introduction

## Summary

## Quick Start
Setup and run the pipeline in 4 steps.

### 1. Prerequisites
#### Hardware
You'll need at least 16GB of available memory to run the pipeline (required for indexing the miRBase reference).

#### Software
Install
[`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)
(>=21.10.3)

Install [`Singualrity`](https://www.sylabs.io/guides/3.0/user-guide/) (see this
[tutorial](https://singularity-tutorial.github.io/01-installation/) )

### 2. Download HEImiRA Pipeline
Download and extract the [latest release](https://github.com/PlantandFoodResearch/HEImiRA-pipeline/releases)
of the pipeline from GitHub.

e.g. on the Linux command line, you can use `curl` and `tar`.
```
curl -L https://github.com/PlantandFoodResearch/HEImiRA-pipeline/releases/download/v1.0/HEImiRA-pipeline-v1.0.tar.gz | tar xzv

cd heimira-pipeline-v1.0
```

### 3. Build BioPandas
BioPandas is a singularity container with BioPython and Pandas.

The build script will build the `biopandas.sif` singularity image inside the
`./singularity` directory.

```
./singularity/build.sh
```

### 4. Run the pipeline
Run the pipeline, with a Nextflow report.

`--input_files` defines the Fastq files to use as input to the pipeline.

`--host_organism` and `--target_organism` can be a *Genus species*, e.g. 'homo
sapiens', or a miRBase organism code, e.g. 'hsa'.  They need to be present in
[miRBase](https://www.mirbase.org).

These parameters (and many more) are all set in the `nextflow.config` file, and
can be overridden at run-time (but make sure to 'single-quote' them), e.g.:

```
nextflow run -with-report heimira-run-report.html main.nf --input_files './my-fastq-data/*.fastq.gz' --host_organism 'homo sapiens' --target_organism 'Zea mays'
```

### 6. Results
Outputs are written to `./results` by default (you can override the `outdir` parameter or change it in the `nextflow.config`).

## Running on HPC (Slurm)
You can run the pipeline on a Slurm cluster with the preset `slurm` Nextflow profile.
You can alter this profile in the `nextflow.config` too.

```
nextflow run -with-report heimira-big-run-report.html -profile slurm main.nf --input_files './my-big-fastq-data/*.fastq.gz'
```

## Testing
You can run the pipeline with the test data as shown below. The `./test-data`
is tiny so the test should run quickly. 

```
nextflow run -with-report heimira-test-report.html main.nf --outdir 'test-results' --input_files './test-data/*.fastq.gz' --host_organism 'homo sapiens' --target_organism 'Malus domestica'
```
