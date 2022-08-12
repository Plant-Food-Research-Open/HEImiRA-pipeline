process clip_3p_adatper {
    tag "clip 3p adatper: $fqfile.SimpleName"
    container = 'https://depot.galaxyproject.org/singularity/cutadapt:2.10--py37hf01694f_1'
    cpus params.max_cpus

    // single-end fastq
    input:
    path fqfile

    output:
    path '*.clip.fq.gz'

    """
    cutadapt -f fastq \
        --trimmed-only \
        -j $task.cpus \
        -a '$params.clip.adapter3p' \
        -m $params.clip.min_read_len \
        -o ${fqfile.SimpleName}.clip.fq.gz \
        $fqfile
    """
}
