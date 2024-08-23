#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    process_read_data(params.input_path)
    process_qc_preprocessing(params.output_name, params.output_path)
    marker_gene_testing(params.output_name, params.output_path, params.MarkerGene)
}

process process_read_data {
    input:
    val input_path

    script:
    """
    python "${params.ScriptDir}read_data.py" --input_path $input_path
    mkdir -p "${params.output_path}"
    mv * "${params.output_path}"
    """
}

process process_qc_preprocessing {
    input:
    val output_name
    val output_path

    script:
    """
    python "${params.ScriptDir}qc_normalization_clustering.py" --output_name $output_name --input_path $output_path
    mv * "${params.output_path}"
    """
}

process marker_gene_testing {
    input:
    val output_name
    val output_path
    val MarkerGene

    script:
    """
    python "${params.ScriptDir}marker_gene_testing.py" --output_name $output_name --input_path $output_path --markergenelist $MarkerGene
    mv * "${params.output_path}"
    """
}
