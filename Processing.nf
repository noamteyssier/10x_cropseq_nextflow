#!/usr/bin/env nextflow

include { KallistoBUS } from "./processes/kallisto_bus.nf"
include { BUStoolsSort } from "./processes/bustools_sort.nf"
include { BUStoolsInspect } from "./processes/bustools_inspect.nf"
include { BUStoolsCount } from "./processes/bustools_count.nf"
include { BUStoolsCorrect } from "./processes/bustools_correct.nf"
include { BuildH5AD } from "./processes/build_h5ad.nf"
include { GeomuxAssign } from "./processes/geomux_assign.nf"
include { MergeAssignments } from "./processes/merge_assignments.nf"

workflow {

    // prepare transcriptome files
    tx_t2g = file(params.kallisto.tx.t2g, checkIfExists: true)
    tx_index = file(params.kallisto.tx.index, checkIfExists: true)
    reads_10x = Channel
        .fromFilePairs (params.data.reads_10x, checkIfExists: true)

    // prepare pcr files
    pcr_t2g = file(params.kallisto.pcr.t2g, checkIfExists: true)
    pcr_index = file(params.kallisto.pcr.index, checkIfExists: true)
    reads_pcr = Channel
        .fromFilePairs (params.data.reads_pcr, checkIfExists: true)

    // prepare whitelist
    whitelist = file(params.whitelist, checkIfExists: true)
    
    // Launch Processes
    ProcessGroups (
        tx_t2g,
        tx_index,
        params.kallisto.tx.tech,
        reads_10x,
        pcr_t2g,
        pcr_index,
        params.kallisto.pcr.tech,
        reads_pcr,
        whitelist,
    )
}

workflow ProcessGroups {
    take:
        tx_t2g
        tx_index
        tx_tech
        reads_10x
        pcr_t2g
        pcr_index
        pcr_tech
        reads_pcr
        whitelist
    main:
        proc_10x = Process10X(
            tx_t2g,
            tx_index, 
            tx_tech, 
            reads_10x,
            whitelist,
        )
        proc_pcr = ProcessPCR(
            pcr_t2g,
            pcr_index,
            pcr_tech,
            reads_pcr,
            whitelist,
        )

        h5ad_ch = proc_10x.h5ad_ch
            .map { it -> [ it[0].replace("_10X", ""), it[1] ] }
        assignment_ch = proc_pcr.assignments_ch
            .map { it -> [ it[0].replace("_PCR", ""), it[1] ] }
        merged_ch = h5ad_ch.join(assignment_ch)
        IntersectAssignments(merged_ch)

}

workflow Process10X {
    take:
        tx_t2g
        tx_index
        tx_tech
        reads_10x
        whitelist
    main:
        bus_10x = KallistoBUS (
            reads_10x,
            tx_index,
            tx_tech,
            "10x",
        )
        correct_10x = BUStoolsCorrect(
            bus_10x.bus_ch,
            whitelist,
            "10x",
        )
        sort_10x = BUStoolsSort(
            correct_10x.correct_bus_ch,
            "10x",
        )
        BUStoolsInspect(
            sort_10x.sorted_bus_ch,
            bus_10x.ec_ch,
            "10x",
        )
        counts_10x = BUStoolsCount(
            sort_10x.sorted_bus_ch,
            bus_10x.transcripts_ch,
            bus_10x.ec_ch,
            tx_t2g,
            "10x",
        )
        h5ad_10x = BuildH5AD(
            counts_10x.mtx_ch,
            counts_10x.barcodes_ch,
            counts_10x.genes_ch,
            "10x",
            "True",
        )
        h5ad_ch = h5ad_10x.h5ad_ch
    emit:
        h5ad_ch
}

workflow ProcessPCR {
    take:
        pcr_t2g
        pcr_index
        pcr_tech
        reads_pcr
        whitelist
    main:
        bus_pcr = KallistoBUS (
            reads_pcr,
            pcr_index,
            pcr_tech,
            "pcr",
        )
        correct_pcr = BUStoolsCorrect(
            bus_pcr.bus_ch,
            whitelist,
            "pcr",
        )
        sort_pcr = BUStoolsSort(
            correct_pcr.correct_bus_ch,
            "pcr",
        )
        BUStoolsInspect(
            sort_pcr.sorted_bus_ch,
            bus_pcr.ec_ch,
            "pcr",
        )
        counts_pcr = BUStoolsCount(
            sort_pcr.sorted_bus_ch,
            bus_pcr.transcripts_ch,
            bus_pcr.ec_ch,
            pcr_t2g,
            "pcr",
        )
        h5ad_pcr = BuildH5AD(
            counts_pcr.mtx_ch,
            counts_pcr.barcodes_ch,
            counts_pcr.genes_ch,
            "pcr",
            "False",
        )
        assignments = GeomuxAssign(
            h5ad_pcr.h5ad_ch,
        )
        assignments_ch = assignments.assignments_ch
    emit:
        assignments_ch
}

workflow IntersectAssignments {
    take:
        merged_ch
    main:
        MergeAssignments(
            merged_ch,
        )
}
