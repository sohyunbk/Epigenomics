#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myProcess {
    script:
    """
    python /home/sb14489/Epigenomics/Spatial/Test/ProduceTable.py
    mv * /scratch/sb14489/
    """

}

workflow {
    myProcess()
}
