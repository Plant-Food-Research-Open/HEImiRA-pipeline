///////////////////// FASTQC ///////////////////////
process fastqc {
    tag "FastQC: $fqfile.SimpleName"
    container 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0'
    // publishDir "${params.outdir}/reports", pattern: '*.html', mode: 'copy'

    cpus params.max_cpus

    input:
    path fqfile

    output:
    path '*.html', emit: reports
    path '*.zip', emit: data

    """
    fastqc --format fastq --threads $task.cpus --noextract --nogroup ${fqfile}
    """
}

///////////////////// MULTIQC /////////////////////
multiqc_config = "$moduleDir/assets/multiqc.yaml"
process multiqc {
    tag "Collate QC results"
    container 'https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0'
    publishDir params.outdir, mode: 'copy'

    input:
    path files
    path config

    output:
    path 'QCsummary.html'
    path 'QCsummary_data'

    """
    multiqc -ip -n QCsummary \
        --module fastqc \
        --config $config \
        --title "$params.title" \
        .
    """
}

workflow qc {
    take:
        file_list
    main:
    fastqc(file_list)
    multiqc(fastqc.out.data.collect(), multiqc_config)
}
