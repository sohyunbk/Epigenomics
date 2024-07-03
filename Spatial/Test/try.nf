#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myProcess {
    output:
    path 'output_table', emit: table

    script:
    """
    python /home/sb14489/Epigenomics/Spatial/Test/ProduceTable.py
    mv output_table /scratch/sb14489/
    """
}

workflow {
    myProcess()
}
