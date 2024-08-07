#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_path = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.output_path = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'
params.ScriptDir = '/home/sb14489/Epigenomics/Spatial/'

workflow {
    process_read_data(params.input_path, params.output_path, params.output_name, params.MarkerGene)
    process_qc_preprocessing("${params.output_path}/adata.h5ad", params.output_path, params.output_name)
    process_normalization("${params.output_path}/adata_qc.h5ad", params.output_path, params.output_name)
    process_clustering("${params.output_path}/adata_norm.h5ad", params.output_path, params.output_name)
    process_marker_gene_testing("${params.output_path}/adata_clustered.h5ad", params.MarkerGene, params.output_path, params.output_name)
}

process process_read_data {
    input:
    path input_path
    path output_path
    val output_name
    path MarkerGene

    output:
    path "${output_path}/adata.h5ad"

    script:
    """
    python "${params.ScriptDir}read_data.py" --input_path $input_path --output_path $output_path --output_name $output_name --MarkerGene $MarkerGene
    """
}

process process_qc_preprocessing {
    input:
    path "${params.output_path}/adata.h5ad"
    path output_path
    val output_name

    output:
    path "${output_path}/adata_qc.h5ad"

    script:
    """
    python "${params.ScriptDir}qc_preprocessing.py" --output_path $output_path --output_name $output_name
    """
}

process process_normalization {
    input:
    path "${params.output_path}/adata_qc.h5ad"
    path output_path
    val output_name

    output:
    path "${output_path}/adata_norm.h5ad"

    script:
    """
    python "${params.ScriptDir}normalization.py" --input_file ${params.output_path}/adata_qc.h5ad --output_file ${params.output_path}/adata_norm.h5ad --output_name $output_name
    """
}

process process_clustering {
    input:
    path "${params.output_path}/adata_norm.h5ad"
    path output_path
    val output_name

    output:
    path "${output_path}/adata_clustered.h5ad"

    script:
    """
    python "${params.ScriptDir}clustering.py" --input_file ${params.output_path}/adata_norm.h5ad --output_file ${params.output_path}/adata_clustered.h5ad
    """
}

process process_marker_gene_testing {
    input:
    path "${params.output_path}/adata_clustered.h5ad"
    path MarkerGene
    path output_path
    val output_name

    script:
    """
    python "${params.ScriptDir}marker_gene_testing.py" --MarkerGene $MarkerGene --output_path $output_path --output_name $output_name --input_file ${params.output_path}/adata_clustered.h5ad
    """
}
