#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process process_read_data {
    input:
    val input_path

    output:
    path("${params.output_path}/read_data_output")

    script:
    """
    python "${params.ScriptDir}read_data.py" --input_path $input_path
    mkdir -p "${params.output_path}/read_data_output"
    mv * "${params.output_path}/read_data_output"
    """
}

process process_qc_preprocessing {
    input:
    path read_data_output

    output:
    path("${params.output_path}/qc_output")

    script:
    """
    python "${params.ScriptDir}qc_normalization_clustering.py" --input_path ${read_data_output}
    mkdir -p "${params.output_path}/qc_output"
    mv * "${params.output_path}/qc_output"
    """
}

process marker_gene_testing {
    input:
    path qc_output
    val MarkerGene

    output:
    path("${params.output_path}/marker_output")

    script:
    """
    python "${params.ScriptDir}marker_gene_testing.py" --input_path ${qc_output} --markergenelist $MarkerGene
    mkdir -p "${params.output_path}/marker_output"
    mv * "${params.output_path}/marker_output"
    """
}

workflow {
    input_path = params.input_path
    output_path = params.output_path
    MarkerGene = params.MarkerGene

    read_data_output = process_read_data(input_path)
    qc_output = process_qc_preprocessing(read_data_output)
    marker_output = marker_gene_testing(qc_output, MarkerGene)
}
