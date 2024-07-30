#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process process_read_data {
    input:
    val input_path

    output:
    path "read_data_output"

    script:
    """
    python "${params.ScriptDir}/read_data.py" --input_path $input_path --output_path read_data_output
    """
}

process process_qc_preprocessing {
    input:
    path "read_data_output"

    output:
    path "qc_output"

    script:
    """
    python "${params.ScriptDir}/qc_normalization_clustering.py" --input_path read_data_output --output_path qc_output
    """
}

process marker_gene_testing {
    input:
    path "qc_output"
    val MarkerGene

    output:
    path "marker_output"

    script:
    """
    python "${params.ScriptDir}/marker_gene_testing.py" --input_path qc_output --markergenelist $MarkerGene --output_path marker_output
    """
}

workflow {
    input_path = params.input_path
    output_path = params.output_path
    MarkerGene = params.MarkerGene

    // Create channels
    read_data_input = Channel.value(input_path)
    read_data_output_channel = Channel.fromPath("read_data_output")
    qc_output_channel = Channel.fromPath("qc_output")
    marker_output_channel = Channel.fromPath("marker_output")

    // Define the workflow steps
    process_read_data(read_data_input).set { read_data_output_channel }
    process_qc_preprocessing(read_data_output_channel).set { qc_output_channel }
    marker_gene_testing(qc_output_channel, MarkerGene).set { marker_output_channel }

    // Move final output to the specified output path
    marker_output_channel.view { path ->
        exec """
        mkdir -p ${output_path}
        mv ${path}/* ${output_path}
        """
    }
}
