process BUStoolsInspect {
    
    publishDir "${params.outdir}/${group}/kb/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(sorted_bus)
    path(ec)
    val(group)

    output:
    path("inspect.json")

    script:
    """
    bustools inspect \
        -o inspect.json \
        -e ${ec} \
        ${sorted_bus}
    """
}
