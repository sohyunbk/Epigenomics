#!/bin/bash
#SBATCH --job-name=scRNA-seq        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=2             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/3_Seurat.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/3_Seurat.%j.err    # Standard error log
#SBATCH --array=0-3                  # Array range

Names=(WTRe1andRe2 WTRe1andRe2_UMI1000 Bif3Re1andRe2 Bif3Re1andRe2_UMI1000)
Re1Objects=(obj_afterDoubletWT_Re1.rds obj_afterDoubletWT_Re1_Gene1000_UMI1000.rds obj_afterDoubletBif3_Re1.rds obj_afterDoubletBif3_Re1_Gene1000_UMI1000.rds)
Re2Objects=(obj_afterDoubletWT_Re2.rds obj_afterDoubletWT_Re2_Gene1000_UMI1000.rds obj_afterDoubletBif3_Re2.rds obj_afterDoubletBif3_Re2_Gene1000_UMI1000.rds)

source activate /home/sb14489/miniconda3/envs/Spatial

Rscript /home/sb14489/Epigenomics/scRNA-seq/4_Harmony.R  \
--WD /scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/ \
--Name "${Names[SLURM_ARRAY_TASK_ID]}" \
--rds1 /scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/"${Re1Objects[SLURM_ARRAY_TASK_ID]}" \
--rds2 /scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/"${Re2Objects[SLURM_ARRAY_TASK_ID]}"
