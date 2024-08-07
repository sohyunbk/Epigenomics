#!/bin/bash
#SBATCH --job-name=scRNA-seq        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=6             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Harmony.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Harmony.%j.err    # Standard error log
#SBATCH --array=0-3                  # Array range

Names=(WTRe1andRe2  WTRe1andRe2_UMI1000  Bif3Re1andRe2  Bif3Re1andRe2_UMI1000)

source activate /home/sb14489/miniconda3/envs/Spatial

Rscript /home/sb14489/Epigenomics/scRNA-seq/5_Annotation.R   \
--WD /scratch/sb14489/4.scRNAseq/2.snRNA-seq/5.MarkerGene/ \
--Name "${Names[SLURM_ARRAY_TASK_ID]}" \
--HarmonyRDS /scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_"${Names[SLURM_ARRAY_TASK_ID]}".rds \
--Markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/SeletedMarkergeneForDotPlot_RemoveUnknown.txt
