#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myProcess {
    output:
    path '/scratch/sb14489/'

    script:
    """
    python /home/sb14489/Epigenomics/Spatial/Test/ProduceTable.py
    """
}

workflow {
    myProcess()
}
