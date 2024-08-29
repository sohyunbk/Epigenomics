process process_read_data {
    output:
    stdout

    """
    #!/home/sb14489/miniconda3/envs/Spatial/bin/python
    import scanpy as sc
    import squidpy as sq
    import pandas as pd
    import os

    # Create the output directory if it does not exist
    os.makedirs("$params.output_path", exist_ok=True)

    # Read the data
    adata = sq.read.visium("$params.input_path")

    pd.set_option('display.max_columns', None)

    # Make variable names unique and calculate QC metrics
    adata.var_names_make_unique()
    adata.var.set_index('gene_ids', inplace=True)
    adata.var.index.name = 'gene_id'

    adata.var['Mt'] = adata.var_names.str.startswith('chrMt')
    adata.var['PT'] = adata.var_names.str.startswith('chrPt')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['Mt', 'PT'], inplace=True)

    # Write the output file to the specified output path
    output_file = os.path.join("$params.output_path", "adata.h5ad")
    adata.write(output_file)
    """
}


process process_qc_preprocessing {
    input:
    stdin

    output:
    stdout

    """
    #!/home/sb14489/miniconda3/envs/Spatial/bin/python
    import scanpy as sc
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os
    print("Current working directory:", os.getcwd())
    # Create the output directory if it does not exist
    adata = sc.read("$params.output_path"+"/adata.h5ad")

    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1], )

    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2]) ## okay line

    sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3], )

    """
}


workflow {
    process_read_data | process_qc_preprocessing | view
    }
