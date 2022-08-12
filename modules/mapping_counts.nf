process mapping_counts {
    tag "Counting mapped reads for $bam"
    container 'https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0'
    publishDir params.outdir

    input:
    path bam

    output:
    path "${bam.getBaseName()}.counts"

    // count no-missmatches only with: grep 'nM:i:0\$'
    """
    samtools view \
    $params.sam_flags \
    -q $params.sam_qual \
    $bam | \
    grep 'nM:i:0\$' | \
    awk '{print \$3}' | \
    sort | uniq -c | \
    awk '{print \$2, \$1}' \
    > ${bam.getBaseName()}.counts

    """
}
