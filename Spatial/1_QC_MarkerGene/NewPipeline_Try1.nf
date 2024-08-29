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
    process_read_data(params.input_path, params.output_path) | process_qc_preprocessing(params.output_name, params.output_path) | marker_gene_testing(params.output_name, params.output_path, params.MarkerGene)
}
