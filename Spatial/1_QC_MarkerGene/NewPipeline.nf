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
    plt.savefig("$params.output_name"+"_QC_Histogram.pdf") ## Save Figure

    #print("Before filtering: Cell - "+ str(cell - adata.n_obs) +"gene - "+str(adata.n_vars))       # check how many genes X cells
    sc.pp.filter_cells(adata, min_counts=100)
    sc.pp.filter_cells(adata, min_genes=50)

    adata = adata[adata.obs["total_counts_Mt"] < 20].copy()
    adata = adata[adata.obs["total_counts_PT"] < 20].copy()
    #print(f"#cells after MT filter: {adata.n_obs}")
    plt.savefig("$params.output_name"+"_QC_Histogram.pdf")
    ## normalization
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    # Manifold embedding and clustering based on transcriptional similarity
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)
    # Save UMAP plot
    plt.rcParams["figure.figsize"] = (4, 4)
    sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save="_" + "$params.output_name")
    plt.rcParams["figure.figsize"] = (8, 8)
    spatial_coords = adata.obsm['spatial'].astype(float)
    adata.obsm['spatial'] = spatial_coords
    sc.pl.spatial(adata, img_key="hires", color=["clusters","total_counts", "n_genes_by_counts"], wspace=0.4, save="_"+"$params.output_name")
    sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[700, 1000, 0, 600],
    alpha=0.5,
    size=1.3,
    save="Magnify_"+"$params.output_name")
    adata.write("$params.output_path"+"/adata_processed.h5ad")
    """
}

process markergene {
    input:
    stdin

    output:
    stdout

    """
    #!/home/sb14489/miniconda3/envs/Spatial/bin/python
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import math
    import os
    import math
    import matplotlib as mpl


    adata = sc.read("$params.input_path"+"/adata_processed.h5ad")
    df = pd.read_csv("$params.MarkerGene", sep='\t')

    gene_list = df['geneID'].tolist()
    #gene_symbols = dict(zip(df['geneID'], df['geneID']+"_"+df['name']))
    gene_symbols = dict(zip(df['geneID'], df['name']))


    gene_ids_in_adata = adata.var.index.values
    filtered_gene_list = [gene for gene in gene_list if gene in gene_ids_in_adata]

    #filtered_gene_list = filtered_gene_list[1:10]

    num_genes = len(filtered_gene_list)
    ncols = 10  # Increase the number of columns to fit more plots in each row
    nrows = math.ceil(num_genes / ncols)

    # Create the spatial plot with larger figure size
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, nrows * 2.5))  # Adjust figsize for larger plots
    axes = axes.flatten()  # Flatten the axes array for easy iteration

    for i, gene in enumerate(filtered_gene_list):
        sc.pl.spatial(
            adata,
            color=gene,
            ax=axes[i],  # Plot on the specific subplot axis
            show=False  # Do not show immediately, we will save it manually
        )
        # Set custom title using gene_symbols dictionary with smaller font size
        if gene in gene_symbols:
            axes[i].set_title(f'{gene}\n{gene_symbols[gene]}', fontsize=8)  # Smaller font size
        else:
            axes[i].set_title(gene, fontsize=8)  # Use gene ID if symbol not found, with smaller font size
        # Customize the scale bar
        for child in axes[i].get_children():
            if isinstance(child, mpl.collections.PatchCollection):
                for path in child.get_paths():
                    if path.vertices.shape[0] == 5:  # Typical of the scale bar
                        path.vertices *= 0.1  # Scale down the size
                        break

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Adjust layout
    plt.tight_layout()

    # Close all figures to avoid the RuntimeWarning
    plt.savefig("MarkerGeneAll_"+"$params.output_name"+".pdf")
    plt.close('all')
    """
    }


workflow {
    process_read_data | process_qc_preprocessing | markergene | view
    }
