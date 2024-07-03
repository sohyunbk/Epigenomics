#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = '/scratch/sb14489/9.spatialRNAseq/1.Rawdata/zma_ear_bif_A'
params.Dir_output = '/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/Bif3_A/'
params.output_name = 'Bif3_A'
params.MarkerGene = '/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker_SelectedMarkerforDotPlot.txt'
params.ScriptDir = '/home/sb14489/Epigenomics/Spatial/'

workflow {
    process_read_data(params.input)
    process_qc_preprocessing(params.output_name, params.Dir_output)
    marker_gene_testing(params.output_name, params.Dir_output, params.MarkerGene)
}

process process_read_data {
    input:
    val input

    script:
    """
    python "${params.ScriptDir}"read_data.py --input_path $input
    mkdir -p "${params.Dir_output}"
    mv * "${params.Dir_output}"
    """
}

process process_qc_preprocessing {
    input:
    val output_name
    val Dir_output
    script:
    """
    python "${params.ScriptDir}qc_normalization_clustering.py" --output_name $output_name --input_path $Dir_output
    mv * "${params.Dir_output}"
    """
}


process marker_gene_testing {
    input:
    val output_name
    val Dir_output
    val MarkerGene
    script:
    """
    python "${params.ScriptDir}qc_normalization_clustering.py" --output_name $output_name --input_path $Dir_output --markergenelist MarkerGene
    mv * "${params.Dir_output}"
    """
}
