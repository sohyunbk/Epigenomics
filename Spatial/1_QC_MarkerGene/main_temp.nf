#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_path = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.output_path = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'
params.ScriptDir = '/home/sb14489/Epigenomics/Spatial/'

workflow {
    process_read_data(params.input_path, params.output_path, params.output_name, params.MarkerGene)
}

process process_read_data {
    input:
    val input_path
    val output_path
    val output_name
    val MarkerGene

    output:
    path "${output_path}/adata.h5ad"

    script:
    """
    python "${params.ScriptDir}read_data.py" --input_path $input_path --output_path $output_path --output_name $output_name --MarkerGene $MarkerGene
    """
}
