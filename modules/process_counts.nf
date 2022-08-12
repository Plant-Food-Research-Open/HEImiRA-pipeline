////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// PROCESS COUNTS ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
process process_counts {
    tag "process alignment counts"
    container 'singularity/biopandas.sif'
    publishDir params.outdir, mode: 'copy', overwrite: true
    cpus params.max_cpus
    // echo true

    input:
    path indexed_bams
    path metadata

    output:
    path "${output_prefix}*.csv"

    shell:
    clean_names = params.clean_sample_names
    output_prefix = params.output_prefix
    template 'process_counts.py'

}

