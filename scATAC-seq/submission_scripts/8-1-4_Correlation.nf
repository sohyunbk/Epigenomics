#!/usr/bin/env nextflow

params.outdir = '/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/1.Correlation'
params.logdir = '/scratch/sb14489/0.log'

// Define standard SLURM options for all processes
executor {
    name = 'slurm'
    queueSize = 2
}

// Process to run Correlation analysis
process CorrelationAnalysis {
    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val sample_id from 0..1

    output:
    path "*.out"
    path "*.err"

    script:
    """
    source activate r_env
    Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Correlation/Correlation_Intergenic2000MostVariationACR_500pbCommonACR.R \\
      --S1_Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${params.S1SampleSparse[sample_id]}" \\
      --S2_Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${params.S2SampleSparse[sample_id]}" \\
      --S1Name A619 \\
      --S2Name "${params.S2SampleName[sample_id]}" \\
      --S1_Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt \\
      --S2_Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${params.S2Metas[sample_id]}" \\
      --ClusterColumnName Ann_v4 \\
      --S1and2_500bpPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${params.CommonPeak[sample_id]}" \\
      --S1and2_500bpInterPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${params.CommonPeakIntergenic[sample_id]}" \\
      --OutPath $PWD \\
      --OutFileName "${params.OutfileName[sample_id]}_pearson" \\
      --CellTypeOrder "${params.CTNameOrder[sample_id]}"
    """
}
