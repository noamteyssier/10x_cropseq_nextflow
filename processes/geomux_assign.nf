process GeomuxAssign {
    
    publishDir "${params.outdir}/pcr/assignments", mode: 'symlink'
    conda "${params.conda.env}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}.tsv"), emit: assignments_ch

    script:
    """
    geomux \
        -i ${h5ad} \
        -o ${sample_id}.tsv \
        -u ${params.geomux.min_umi} \
        -c ${params.geomux.min_cells} \
        -j ${params.n_threads}
    """
}
