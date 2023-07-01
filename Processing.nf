#!/usr/bin/env nextflow

include { KallistoBUS } from "./processes/kallisto_bus.nf"
include { BUStoolsSort } from "./processes/bustools_sort.nf"
include { BUStoolsInspect } from "./processes/bustools_inspect.nf"
include { BUStoolsCount } from "./processes/bustools_count.nf"
include { BuildH5AD } from "./processes/build_h5ad.nf"
include { GeomuxAssign } from "./processes/geomux_assign.nf"

workflow {

  pcr_t2g = file(params.kallisto.pcr.t2g, checkIfExists: true)
  pcr_index = file(params.kallisto.pcr.index, checkIfExists: true)
  tx_t2g = file(params.kallisto.tx.t2g, checkIfExists: true)
  tx_index = file(params.kallisto.tx.index, checkIfExists: true)

  reads_pcr = Channel
    .fromFilePairs (params.data.reads_pcr, checkIfExists: true)

  // reads_10x = Channel
  //   .fromFilePairs (params.data.reads_10x, checkIfExists: true)

  // PCR Pipeline
  pcr_bus = KallistoBUS (
    reads_pcr,
    pcr_index,
    params.kallisto.pcr.tech,
    "pcr",
  )
  pcr_sort = BUStoolsSort(
    pcr_bus.bus_ch,
    "pcr",
  )
  BUStoolsInspect(
    pcr_sort.sorted_bus_ch,
    pcr_bus.ec_ch,
    "pcr",
  )
  pcr_counts = BUStoolsCount(
    pcr_sort.sorted_bus_ch,
    pcr_bus.transcripts_ch,
    pcr_bus.ec_ch,
    pcr_t2g,
    "pcr",
  )
  pcr_h5ad = BuildH5AD(
    pcr_counts.mtx_ch,
    pcr_counts.barcodes_ch,
    pcr_counts.genes_ch,
    "pcr",
    "False",
  )
  pcr_assignments = GeomuxAssign(
    pcr_h5ad.h5ad_ch,
  )
}



// process BuildH5AD {
//     
//     publishDir "${params.outdir}/counts/${sample_id}", mode: 'symlink'
//     conda "${params.scanpy_env}"

//     input:
//     tuple val(sample_id), path(mtx)
//     path(barcodes)
//     path(genes)

//     output:
//     tuple val(sample_id), path("${sample_id}.h5ad"), emit: h5ad_ch

//     script:
//     """
//     #!/usr/bin/env python3

//     import scanpy as sc
//     import anndata as ad

//     # Load data
//     adata = sc.read_mtx("${mtx}")
//     barcodes = open("${barcodes}").read().splitlines()
//     genes = open("${genes}").read().splitlines()

//     # Set index
//     adata.obs.index = barcodes
//     adata.var.index = genes

//     # Calculate QC metrics
//     sc.pp.calculate_qc_metrics(adata, inplace=True)

//     # Add SampleID to obs
//     adata.obs["sample_id"] = "${sample_id}"

//     # Write to file
//     adata.write_h5ad("${sample_id}.h5ad")
//     """
// }

// process PlotQC {
//     
//     publishDir "${params.outdir}/qc", mode: 'symlink'
//     conda "${params.scanpy_env}"

//     input:
//     tuple val(sample_id), path(h5ad)

//     output:
//     path("${sample_id}_rankplot.svg")

//     script:
//     """
//     #!/usr/bin/env python3

//     import numpy as np
//     import scanpy as sc
//     import matplotlib.pyplot as plt

//     # Load data
//     adata = sc.read_h5ad("${h5ad}")

//     # Generate Axis
//     knee = np.sort(np.asarray(adata.X.sum(axis=1)).ravel())[::-1]
//     ranks = np.arange(knee.size)

//     # Fraction of UMIs kept with filter
//     total_cells = adata.shape[0]
//     kept_cells = (knee > ${params.cell_umi_threshold}).sum()
//     frac = kept_cells / total_cells

//     # Plot
//     fig, ax = plt.subplots(figsize=(10, 7))

//     ax.loglog(knee, ranks, linewidth=5, color="g")
//     ax.axvline(${params.cell_umi_threshold}, color="r", linestyle="--", linewidth=3)

//     ax.set_xlabel("UMI Counts")
//     ax.set_ylabel("Set of Barcodes")
//     ax.set_title("${sample_id} - Fraction of Droplets Kept: {kept} / {total} ({frac:.4f})".format(total=total_cells, kept=kept_cells, frac=frac))

//     plt.grid(True, which="both")
//     plt.tight_layout()
//     plt.savefig("${sample_id}_rankplot.svg")
//     """
// }

// process FilterNormLog {
//     
//     publishDir "${params.outdir}/counts/${sample_id}", mode: 'symlink'
//     conda "${params.scanpy_env}"

//     input:
//     tuple val(sample_id), path(h5ad)

//     output:
//     path("${sample_id}.filt.h5ad"), emit: h5ad_filt_ch
//     val(sample_id), emit: sample_id_ch

//     script:
//     """
//     #!/usr/bin/env python3

//     import scanpy as sc

//     # Load data
//     adata = sc.read_h5ad("${h5ad}")

//     # Filter cells
//     sc.pp.filter_cells(adata, min_counts=${params.cell_umi_threshold})

//     # Filter genes
//     sc.pp.filter_genes(adata, min_counts=${params.gene_cell_threshold})

//     # Set Raw Data
//     adata.raw = adata

//     # Normalize and Log
//     sc.pp.normalize_total(adata, target_sum=1e4)
//     sc.pp.log1p(adata)

//     # Write to file
//     adata.write_h5ad("${sample_id}.filt.h5ad")
//     """
// }

// process Merge {

//     publishDir "${params.outdir}/counts", mode: 'symlink'
//     conda "${params.scanpy_env}"

//     input:
//     path(adata_filt_ch)

//     output:
//     path("adata.h5ad"), emit: h5ad_ch

//     script:
//     """
//     #!/usr/bin/env python3

//     import scanpy as sc

//     def build_g2s(t2g_path: str) -> dict:
//         used = dict()
//         g2s = dict()
//         for line in open(t2g_path):
//             record = line.strip().split('\t')
//             tx = record[1]
//             sym = record[2]
//             if tx not in g2s:
//                 if sym == "":
//                     g2s[tx] = tx
//                 else:
//                     if sym not in used:
//                         used[sym] = 0
//                     else:
//                         used[sym] += 1
//                         sym += f"_{used[sym]}"
//                     g2s[tx] = sym
//         return g2s

//     def concatenate_h5ads(path_list: list):
//         adata = sc.concat(
//             [sc.read_h5ad(i) for i in path_list],
//             join='outer',
//         )
//         adata.obs_names_make_unique()
//         return adata


//     path_list = [i for i in "${adata_filt_ch}".split(" ")]
//     adata = concatenate_h5ads(path_list)
//     g2s = build_g2s("${params.t2g}")
//     
//     adata.var["ensembl_id"] = adata.var.index
//     adata.var["gene_symbol"] = adata.var["ensembl_id"].map(g2s)
//     adata.var.index = adata.var["gene_symbol"]
//     adata.var_names_make_unique()

//     adata.write_h5ad("adata.h5ad")
//     """
// }

// process SingleCell {

//     publishDir "${params.outdir}/analysis", mode: 'symlink'
//     conda "${params.scanpy_env}"

//     input:
//     path(merged_h5ad)

//     output:
//     path(merged_h5ad), emit: sc_h5ad_ch

//     script:
//     """
//     #!/usr/bin/env python3

//     import scanpy as sc

//     # Load data
//     adata = sc.read_h5ad("${merged_h5ad}")

//     # Calculate Highly Variable Genes
//     sc.pp.highly_variable_genes(adata, n_top_genes=2000)

//     # Perform PCA
//     sc.tl.pca(adata)

//     # Perform Nearest Neighbor Graph
//     sc.pp.neighbors(adata)

//     # Perform UMAP
//     sc.tl.umap(adata)

//     # Perform Leiden Clustering
//     sc.tl.leiden(adata)

//     # Write to file
//     adata.write_h5ad("${merged_h5ad}")
//     """

// }
