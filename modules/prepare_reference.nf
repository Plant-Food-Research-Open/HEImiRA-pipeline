///////////////////////////////////////////////////////////////////////////////
////////////////////////////// PREPARE REFERENCE //////////////////////////////
///////////////////////////////////////////////////////////////////////////////
process make_reference {
    tag "Make collapsed reference"
    container 'singularity/biopandas.sif'
    publishDir params.outdir, mode: 'copy'

    input:
    path ref
    path taxa

    output:
    path out_fasta, emit: fasta
    path out_table, emit: data

    shell:
    out_fasta = params.reference.out.fasta
    out_table = params.reference.out.table
    host_organism = params.host_organism
    target_organism = params.target_organism

    template 'prepare_reference.py'

}

workflow prepare_reference {
    main:
    make_reference(
        params.reference.seqs,
        params.reference.taxa
    )

    emit:
    fasta = make_reference.out.fasta
    data = make_reference.out.data
}
