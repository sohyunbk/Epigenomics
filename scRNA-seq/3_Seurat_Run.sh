#!/bin/bash
#SBATCH --job-name=scRNA-seq        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=8             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/3_Seurat.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/3_Seurat.%j.err    # Standard error log
#SBATCH --array=0-7                  # Array range

Names=(WT_Re1_Gene1000_UMI1000 WT_Re2_Gene1000_UMI1000 Bif3_Re1_Gene1000_UMI1000 Bif3_Re2_Gene1000_UMI1000 WT_Re1 WT_Re2 Bif3_Re1 Bif3_Re2)
UMICuts=(1000 1000 1000 1000 0 0 0 0)
GeneCuts=(1000 1000 1000 1000 0 0 0 0)
InputFiles=(Sohyun-wt-1 Sohyun-wt-2 Sohyun-bif3-1 Sohyun-bif3-2 Sohyun-wt-1 Sohyun-wt-2 Sohyun-bif3-1 Sohyun-bif3-2)

source activate /home/sb14489/miniconda3/envs/Spatial

Rscript /home/sb14489/Epigenomics/scRNA-seq/3_Seurat_QC.R \
--WD /scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/ \
--Name "${Names[SLURM_ARRAY_TASK_ID]}" \
--UMICut "${UMICuts[SLURM_ARRAY_TASK_ID]}" \
--GeneCut "${GeneCuts[SLURM_ARRAY_TASK_ID]}" \
--InputDir /scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/"${InputFiles[SLURM_ARRAY_TASK_ID]}"/outs/raw_feature_bc_matrix/
