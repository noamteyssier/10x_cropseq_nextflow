process BUStoolsCount {

    publishDir "${params.outdir}/${group}/counts/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(sorted_bus)
    path(transcripts)
    path(t2g)
    path(ec)
    val(group)

    output:
    tuple val(sample_id), path("data.mtx"), emit: mtx_ch
    path "data.barcodes.txt", emit: barcodes_ch
    path "data.genes.txt", emit: genes_ch

    script:
    """
    bustools count \
        -o data \
        -t ${transcripts} \
        -g ${t2g} \
        -e ${ec} \
        --genecounts \
        ${sorted_bus}
    """
}
