#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myProcess {
    output:
    path 'output', emit: my_output

    script:
    """
    python ProduceTable.py
    """
}

workflow {
    myProcess()
}
