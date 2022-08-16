# HEImiRA Pipeline
_Host environment intersection miRNA annotation_ pipeline

## Introduction
HEImiRA (Host Environment Intersection miRNA Annotation) pipeline is built with
Nextflow using singularity containers from the [Galaxy Project] and a custom
python container, to annotate miRNAs represented in sRNA-Seq data against the
[miRBase](https://www.mirbase.org/) miRNA mature reference sequences.

## Summary
### Input
sRNA-Seq reads, `FASTQ` format, optionally gzipped.

### Preprocessing
Preparing the sRNA-Seq read data for alignment.

 - Pre and post quality checking - [FASTQC + MultiQC](modules/qc.nf)
 - Adaptor clipping - [cutadapt](modules/clip.nf)
 - Filtering - [BBDUK](modules/filter.nf)

### miRBase Tagging
Prepare the miRBase reference for alignment - [`prepare_reference.py`](templates/prepare_reference.py) (HEImiRA script).

 - miRBase reference sequences collapsed
 - Taxonomic groups ascribed to collapsed sequences
 - Sequences tagged by `target`, `host`, `environment`
 - Create collapsed-sequence reference FASTA `heimeria_collapsed_reference.fa.gz`
 - Create collapsed-sequence + taxonomy table `heimeria_metadata.csv.gz`

### Alignment
Align preprocessed sRNA-Seq reads to the HEImiRA collapsed miRBase reference.

 - Align the reads to the reference - [STAR + Samtools](modules/map_reads.nf)

### Analysis
Apply HEImiRA tags and combine taxonomic metadata with alignment counts - [`process_counts.py`](templates/process_counts.py) (HEImiRA script).

 - Generate alignment counts
 - Apply HEImiRA `target`, `host`, and `environment` tags
 - Combine taxonomic metadata

### Output
HEImiRA pipeline outputs are in CSV format compatible with [Python Pandas](https://pandas.pydata.org/)
for easy direct use or downstream analysis.

The 3 output tables contain the HEImiRA tags and taxonomic metadata for each
reference sequence. The alignment counts are represented differently in each
table, as follows.

 - **HEImiRA counts summary table**
   - raw counts
 - **HEImiRA normalised counts summary table**
   - counts normalised by total counts per sample
 - **HEImiRA host-normalised counts summary table**
   - counts normalised by total host-tagged counts per sample

## Quick Start
Setup and run the pipeline in 4 steps.

### 1. Prerequisites
#### Hardware
You'll need at least 16GB of available memory to run the pipeline (required for
indexing the miRBase reference).

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

### Results
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
