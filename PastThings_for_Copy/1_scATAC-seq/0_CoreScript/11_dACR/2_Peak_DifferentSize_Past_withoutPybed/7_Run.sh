#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11_Fisher.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11_Fisher.%j.err    # Standard error log
#SBATCH --array=0-12

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)
ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/r_env

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/11_dACR/2_Peak_DifferentSize/7_FisherExactTest.R "${ClusterN[SLURM_ARRAY_TASK_ID]}"
