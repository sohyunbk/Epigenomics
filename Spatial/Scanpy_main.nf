#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_path = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.output_path = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'

workflow {
    adata = process_read_data(params.input_path, params.output_path, params.output_name, params.MarkerGene)
    adata_qc = process_qc_preprocessing(adata, params.output_path, params.output_name)
    adata_norm = process_normalization(adata_qc)
    adata_clustered = process_clustering(adata_norm)
    process_marker_gene_testing(adata_clustered, params.MarkerGene, params.output_path, params.output_name,adata_clustered)
}

process process_read_data {
    input:
    path input_path
    path output_path
    path output_name
    path MarkerGene

    output:
    path "adata.h5ad" into adata

    script:
    """
    python read_data.py --input_path $input_path --output_path $output_path --output_name $output_name --MarkerGene $MarkerGene
    """
}

process process_qc_preprocessing {
    input:
    path output_path
    path output_name

    output:
    path "adata_qc.h5ad" into adata_qc

    script:
    """
    python qc_preprocessing.py --output_path $output_path --output_name $output_name
    """
}

process process_normalization {
    input:
    path "adata_qc.h5ad"

    output:
    path "adata_norm.h5ad" into adata_norm

    script:
    """
    python normalization.py --input_file adata_qc.h5ad
    """
}

process process_clustering {
    input:
    path "adata_norm.h5ad"

    output:
    path "adata_clustered.h5ad" into adata_clustered

    script:
    """
    python clustering.py --input_file adata_norm.h5ad
    """
}

process process_marker_gene_testing {
    input:
    path MarkerGene
    path output_path
    path output_name
    path "adata_clustered.h5ad"
    script:
    """
    python marker_gene_testing.py --MarkerGene $MarkerGene --output_path $output_path --output_name $output_name --input_file adata_clustered.h5ad
    """
}
