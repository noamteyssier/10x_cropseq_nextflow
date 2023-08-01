process MergeAssignments {
    
    publishDir "${params.outdir}/adata", mode: 'copy'
    conda "${params.conda.env}"

    input:
    tuple val(sample_id), path(adata), path(assignment)

    output:
    path("${sample_id}.h5ad")

    script:
    """
    #!/usr/bin/env python3
    
    import scanpy as sc
    import pandas as pd

    # load adata
    adata = sc.read_h5ad(
        "${adata}"
    )

    # load assignment
    assignment = pd.read_csv(
        "${assignment}",
        sep='\t',
    )

    # parse assignment name for singletons
    singles = assignment[assignment.moi == 1].copy()
    singles["guide"] = singles.assignment.apply(
        lambda x: x.replace("['", "").replace("']", "")
    )

    # build cell to guide mapper
    mapper = singles.loc[:, ["cell_id", "guide"]]\
        .set_index("cell_id")\
        .to_dict()["guide"]

    # assign guides to cells
    adata.obs["guide"] = adata.obs.index.map(mapper)

    # save adata
    adata.write_h5ad("${sample_id}.h5ad")
    """

}
