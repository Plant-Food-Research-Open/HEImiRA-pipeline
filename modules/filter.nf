process bbduk_filter {
    tag "Filtering: $fqfile.SimpleName"
    container = 'https://depot.galaxyproject.org/singularity/bbmap:38.90--he522d1c_1'
    cpus params.max_cpus

    input:
    path fqfile

    output:
    path "${fqfile.SimpleName}.filter.fq.gz", emit: reads
    path "${fqfile.SimpleName}.stats", emit: stats

    """
    bbduk.sh \
    threads=${task.cpus} \
    in=$fqfile \
    ref=adapters \
    k=$params.filter.k \
    maxns=$params.filter.maxns \
    overwrite=true \
    out=${fqfile.SimpleName}.filter.fq.gz \
    stats=${fqfile.SimpleName}.stats
    """
}


