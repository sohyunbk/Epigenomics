#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process runPipeline {
    // Specify the container, if using one
    // container 'python:3.8'
    script:
    """
    python "${params.ScriptDir}/AllPipeline.py" --input_path params.input_path --output_path params.output_path --output_name params.output_name --markergenelist params.MarkerGene
    """
}

workflow {
    runPipeline()
}
