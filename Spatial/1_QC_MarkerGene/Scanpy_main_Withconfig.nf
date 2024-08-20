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
    path "qc_data_output"

    script:
    """
    python "${params.ScriptDir}/qc_normalization_clustering.py" --input_path read_data_output --output_name qc_data_output
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

    // Step 1: Read data
    read_data_out = process_read_data(input_path)

    // Step 2: QC preprocessing
    qc_out = process_qc_preprocessing(read_data_out)

    // Step 3: Marker gene testing
    marker_out = marker_gene_testing(qc_out, MarkerGene)

    // Step 4: Move final output to the specified output path
    marker_out.view { path ->
        exec """
        mkdir -p ${output_path}
        mv ${path}/* ${output_path}
        """
    }
}
