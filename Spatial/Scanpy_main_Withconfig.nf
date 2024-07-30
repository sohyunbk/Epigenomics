#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process process_read_data {
    input:
    val input_path

    output:
    path("read_data_output") into read_data_out

    script:
    """
    python "${params.ScriptDir}read_data.py" --input_path $input_path
    mkdir -p read_data_output
    mv * read_data_output
    """
}

process process_qc_preprocessing {
    input:
    path read_data_out

    output:
    path("qc_output") into qc_out

    script:
    """
    python "${params.ScriptDir}qc_normalization_clustering.py" --input_path read_data_output
    mkdir -p qc_output
    mv * qc_output
    """
}

process marker_gene_testing {
    input:
    path qc_out
    val MarkerGene

    output:
    path("marker_output") into marker_out

    script:
    """
    python "${params.ScriptDir}marker_gene_testing.py" --input_path qc_output --markergenelist $MarkerGene
    mkdir -p marker_output
    mv * marker_output
    """
}

workflow {
    input_path = params.input_path
    output_path = params.output_path
    MarkerGene = params.MarkerGene

    read_data_out = process_read_data(input_path)
    qc_out = process_qc_preprocessing(read_data_out)
    marker_out = marker_gene_testing(qc_out, MarkerGene)

    // Move final output to the specified output path
    marker_out.view { path ->
        exec """
        mkdir -p ${output_path}
        mv ${path}/* ${output_path}
        """
    }
}
