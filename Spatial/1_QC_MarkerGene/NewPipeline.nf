process process_read_data {
    input:
    path x
    path y
    script:
    """"
    python "${params.ScriptDir}"read_data.py --input_path $x --output_path $y
    mkdir -p "${params.Dir_output}"
    """
}

process process_qc_preprocessing {
    input:
    val outputname
    path Diroutput
    script:
    """
    python "${params.ScriptDir}qc_normalization_clustering.py" --output_name "${params.output_path}"$outputname --input_path $Diroutput
    mv * "${params.Dir_output}"
    """
}


process marker_gene_testing {
    input:
    val outputname
    path Diroutput
    path Marker
    script:
    """
    python "${params.ScriptDir}marker_gene_testing.py" --output_name "${params.output_path}"$outputname --input_path $Diroutput --markergenelist $Marker
    """
}

workflow {
    // Create channels from parameters
    Channel
        .value(params.input)
        .set { input_ch }

    Channel
        .value(params.output_path)
        .set { output_ch }

    Channel
        .value(params.output_name)
        .set { outputname_ch }

    Channel
        .value(params.MarkerGene)
        .set { marker_ch }

    // Call processes with channels
    process_read_data(input_ch, output_ch)
    process_qc_preprocessing(outputname_ch, output_ch)
    marker_gene_testing(outputname_ch, output_ch, marker_ch)
}
