#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myProcess {
    output:
    path '/scratch/sb14489'

    script:
    """
    python ProduceTable.py
    """
}
