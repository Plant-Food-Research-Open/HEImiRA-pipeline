# HEImiRA Pipeline
_Host environment intersection miRNA annotation_ pipeline

## Introduction

## Summary

## Quick Start

### 1.  Prerequisites
#### Hardware
You'll need at least 16GB of available memory to run the pipeline (required for indexing the miRBase reference).

#### Software
Install
[`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)
(>=21.10.3)

Install [`Singualrity`](https://www.sylabs.io/guides/3.0/user-guide/) (see this
[tutorial](https://singularity-tutorial.github.io/01-installation/) )

### 2. Get HEImiRA Pipeline
Clone it from GitHub.

```
git clone <URL> heimira-pipeline

cd heimira-pipeline
```

### 4. Build BioPandas
BioPandas is a singularity container with BioPython and Pandas.

The build script will build the `biopandas.sif` singularity image inside the
`./singularity` directory.

```
./singularity/build.sh
```

### 5. Run the pipeline
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
