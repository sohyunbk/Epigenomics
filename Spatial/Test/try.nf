#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myProcess {
    output:
    path 'output', emit: my_output

    script:
    """
    python /home/sb14489/Epigenomics/Spatial/Test/ProduceTable.py
    """
}

workflow {
    myProcess()
}
