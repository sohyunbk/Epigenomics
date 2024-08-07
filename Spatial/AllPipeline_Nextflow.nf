#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process runPipeline {
    // Specify the container, if using one
    // container 'python:3.8'

    input:
    file inputPath from params.input_path
    file outputPath from params.output_path
    val outputName from params.output_name
    file markerGene from params.MarkerGene
    val scriptDir from params.ScriptDir

    script:
    """
    python "${scriptDir}/AllPipeline.py" --input_path ${inputPath} --output_path ${outputPath} --output_name ${outputName} --markergenelist ${markerGene}
    """
}

workflow {
    runPipeline()
}
