#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/////////////////// INPUT FILES ///////////////////
input_files = Channel.fromPath(params.input_files)

///////////////////// MODULES /////////////////////
// QC module - Raw
include {qc as rawqc} from './modules/qc' addParams(
    outdir: "${params.outdir}/qc/raw",
    title: "Raw FastQC Results"
    )
// QC module - Post
include {qc as postqc} from './modules/qc' addParams(
    outdir: "${params.outdir}/qc/post",
    title: "Post FastQC Results"
    )
// Adapter clipping module
include {clip_3p_adatper as clip} from './modules/clip'
// Filter module
include {bbduk_filter as filter_reads} from './modules/filter'
// Collapse and index miRBase reference with metadata
include {prepare_reference} from './modules/prepare_reference'
// STAR mapping module
include {map_reads} from './modules/map_reads'
// Process counts, counts and summary
include {process_counts} from './modules/process_counts'

///////////////////// WORKFLOW /////////////////////
workflow {
    // Raw Data QC
    rawqc(input_files)
    // Clip adapters
    clip(input_files)
    // Filter out known contaminants
    filter_reads(clip.out)
    // post QC
    postqc(filter_reads.out.reads)
    // Prepare miRNA reference sequences
    prepare_reference()
    // Map reads to prepared reference
    map_reads(filter_reads.out.reads, prepare_reference.out.fasta)
    // Process counts
    process_counts(
        // Collect all bam, bai tuples
        map_reads.out.indexed_bam.collect(),
        prepare_reference.out.data
    )
}

