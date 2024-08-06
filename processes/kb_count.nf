process KBCount {
    maxForks 3

    publishDir "${params.outdir}/${group}/counts/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)
    path(intron_t2g)
    path(cdna_t2g)
    path(t2g)
    path(index)
    val(tech)

    output:
    tuple val(sample_id), val(tech), path("${sample_id}.h5ad"), emit: h5ad_ch

    script:
    """
    kb count \
        --h5ad \
        -i ${index} \
        -g ${t2g} \
        -x ${tech} \
        --workflow lamanno \
        -c1 ${cdna_t2g} \
        -c2 ${intron_t2g} \
        -o ${sample_id} \
        --filter bustools \
        -t ${params.n_threads} \
        ${reads[0]} ${reads[1]}
    """
}