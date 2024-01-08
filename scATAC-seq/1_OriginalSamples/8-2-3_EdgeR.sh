#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11.dACR.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11.dACR.%j.err    # Standard error log
#SBATCH --array=0-13

## 14 number of cell types
ClusterN=(L1atFloralMeristem  Unknown1 FloralMeristem_SuppressedBract  PhloemPrecursor  Unknown2
G2_M ProcambialMeristem_ProtoXylem_MetaXylem  Unknown_Sclerenchyma IM-OC ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
Unknown_lowFRiP L1 SPM-base_SM-base  XylemParenchyma_PithParenchyma)

ml Anaconda3/2023.09-0
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/dACR/Differential_ACR_PseudoBulk.R \
 --Sparse_S1 /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_PeaksCount_byA619Barcode.txt \
 --Sparse_S2 /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_PeaksCount_byBif3Barcode.txt \
 --Meta_S1 /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt  \
 --Meta_S2 /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV4_metadata.txt \
 --CellType "${ClusterN[SLURM_ARRAY_TASK_ID]}" \
 --IntergenicPeakFile /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_MergedDifferentSizePeak/A619Bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_MergedPeak_Intergenic.bed  \
 --OutputDir /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4
