#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_path = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.output_path = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'

workflow {
    process_read_data()
    process_qc_preprocessing()
    process_normalization()
    process_clustering()
    process_visualization()
    process_marker_gene_testing()
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
