process BUStoolsCorrect {

    publishDir "${params.outdir}/${group}/kb/${sample_id}", mode: "symlink"

    input:
    tuple val(sample_id), path(sorted_bus)
    path(whitelist)
    val(group)

    output:
    tuple val(sample_id), path("output.correct.bus"), emit: correct_bus_ch

    script:
    """
    bustools correct \
        -w ${whitelist} \
        -o output.correct.bus \
        ${sorted_bus}
    """
}