#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.Dir = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.output_path = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'
params.ScriptDir = '/home/sb14489/Epigenomics/Spatial/'

workflow {
    process_read_data(params.Dir)
    process_qc_preprocessing(params.output_name, params.Dir)
}

process process_read_data {
    input:
    val Dir

    script:
    """
    python "${params.ScriptDir}"read_data.py --input_path $Dir
    mkdir -p "${params.output_path}"
    mv * "${params.output_path}"
    """
}

process process_qc_preprocessing {
    input:
    val output_name
    val Dir
    script:
    """
    python "${params.ScriptDir}qc_preprocessing.py" --output_name $output_name --input_path $Dir
    """
}
