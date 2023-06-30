process BUStoolsSort {
    
    publishDir "${params.outdir}/${group}/kb/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(unsorted_bus)
    val(group)

    output:
    tuple val(sample_id), path("output.sorted.bus"), emit: sorted_bus_ch

    script:
    """
    bustools sort \
        -t ${params.n_threads} \
        -m ${params.max_mem} \
        -o output.sorted.bus \
        ${unsorted_bus}
    """
}
