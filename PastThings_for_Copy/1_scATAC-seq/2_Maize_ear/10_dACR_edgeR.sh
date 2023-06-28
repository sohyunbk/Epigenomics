#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11_dACR.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11_dACR.%j.err    # Standard error log
#SBATCH --array=0-11

CellType=(CalloseRelated FloralMeristem_SuppressedBract G2_M IM_SPM_SM IM-OC L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)


ml Anaconda3/2020.02
source activate r4-base

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/11_Differential_ACR_SubmitJob.R "${CellType[SLURM_ARRAY_TASK_ID]}"
