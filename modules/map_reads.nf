process star_index {
    tag "STAR: indexing ${ref.getName()}"
    container 'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0'
    cpus params.max_cpus
    memory params.max_memory

    input:
    path ref
    val refsize

    output:
    path index_dir

    script:
    log2refsize = 31 - Integer.numberOfLeadingZeros( refsize as Integer )
    genomeSAindexNbases = Math.min(14.0, log2refsize/2.0 - 1.0) as Integer

    index_dir = 'STAR_index'
    """
    mkdir $index_dir && \
    STAR --runMode genomeGenerate \
    --runThreadN $task.cpus \
    --limitGenomeGenerateRAM ${task.memory.toBytes()} \
    --genomeSAindexNbases $genomeSAindexNbases \
    --genomeFastaFiles $ref \
    --genomeDir $index_dir
    """
}

process star_align {
    tag "Mapping $readfq"
    container 'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0'
    cpus params.max_cpus
    memory params.max_memory
    // Override with extra time for STAR alignment
    time '6h'

    input:
    path readfq
    path index

    output:
    path "${outprefix}.bam"

    script:
    outprefix = readfq.getSimpleName()
    // 80% of available mem for bam sorting (bytes)
    bam_mem_limit = task.memory.toBytes() *0.8 as Integer

    """
    STAR --runMode alignReads \
    --runThreadN $task.cpus \
    --readFilesCommand zcat \
    --alignEndsType EndToEnd \
    --twopassMode None \
    --alignIntronMax 1 \
    --outFilterScoreMin $params.sam_qual \
    --outFilterMismatchNmax $params.alignment_mismatch \
    --outSAMtype BAM SortedByCoordinate \
    --outStd BAM_SortedByCoordinate \
    --limitBAMsortRAM ${bam_mem_limit} \
    --genomeDir $index \
    --readFilesIn $readfq \
    --outFileNamePrefix $outprefix \
    > ${outprefix}.bam
    """
}

process ref_size {
    tag "Getting size of reference ${ref.getName()}"

    input:
    path ref

    output:
    env REFSIZE

    """
    REFSIZE=\$( \
    gzip -dc < $ref | \
    awk 'BEGIN{ ORS=""} /^[^>]/ {print \$0}' | \
    wc -c \
    )

    """
}

process decompress {
    tag "Decompressing $gzfile"

    input:
    path gzfile

    output:
    path "${gzfile.getBaseName()}"

    """
    gzip -dc < $gzfile > ${gzfile.getBaseName()}
    """
}

process samtools_index {
    tag "Indexing bam: $bam"
    container 'https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0'
    cpus params.max_cpus
    // publishDir "${params.outdir}/alignments"

    input:
    path bam

    output:
    tuple(path(bam), path(bai))

    script:
    bai = "${bam}.bai"
    """
    samtools index -@ $task.cpus $bam $bai
    """
}

workflow map_reads {
    take:
    reads
    reference

    main:
    ref_size(reference)
    star_index(decompress(reference), ref_size.out)
    star_align(reads, star_index.out)
    samtools_index(star_align.out)

    // alignments = star_align.out
        // .join(samtools_index.out.map{[ it.getBaseName(), it]})
    // alignments.view()
    emit:
    indexed_bam = samtools_index.out
}
