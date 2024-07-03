#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_path = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.output_path = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'

workflow {
    input_data = Channel.fromPath(params.input_path)
    marker_gene = Channel.fromPath(params.MarkerGene)

    process_read_data(input_data, params.output_path, params.output_name, marker_gene)
    process_qc_preprocessing(params.output_path, params.output_name)
    process_normalization()
    process_clustering()
    process_marker_gene_testing(marker_gene, params.output_path, params.output_name)
}

process process_read_data {
    input:
    path input_path
    path output_path
    path output_name
    path MarkerGene

    output:
    path "adata.h5ad" into adata_channel

    script:
    """
    python read_data.py --input_path $input_path --output_path $output_path --output_name $output_name --MarkerGene $MarkerGene
    """
}

process process_qc_preprocessing {
    input:
    path adata_channel
    path output_path
    path output_name

    output:
    path "adata_qc.h5ad" into adata_qc_channel

    script:
    """
    python qc_preprocessing.py --output_path $output_path --output_name $output_name
    """
}

process process_normalization {
    input:
    path adata_qc_channel

    output:
    path "adata_norm.h5ad" into adata_norm_channel

    script:
    """
    python normalization.py --input_file adata_qc.h5ad
    """
}

process process_clustering {
    input:
    path adata_norm_channel

    output:
    path "adata_clustered.h5ad" into adata_clustered_channel

    script:
    """
    python clustering.py --input_file adata_norm.h5ad
    """
}

process process_marker_gene_testing {
    input:
    path adata_clustered_channel
    path MarkerGene
    path output_path
    path output_name

    script:
    """
    python marker_gene_testing.py --MarkerGene $MarkerGene --output_path $output_path --output_name $output_name --input_file adata_clustered.h5ad
    """
}
