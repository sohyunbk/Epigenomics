#!/bin/bash
#SBATCH --job-name=BulkRe_Corr        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=40:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/BulkRe_Corr.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/BulkRe_Corr.%j.err    # Standard error log
#SBATCH --array=0-2


Re1=(1_A619 2_rel2 3_bif3)
Re2=(1_A619_2 2_rel2_2 3_bif3_2)
OutFileNames=(A619Re1andRe2_Union rel2Re1andRe2 Bif3Re1andRe2_Union)

ml Anaconda3/2022.10
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Viualization/Replicate1and2Correlation_UnionPeaks.R \
  --Re1_BulkPeak /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/"${Re1[SLURM_ARRAY_TASK_ID]}"/macs2_temp/bulk_peaks_peaks.narrowPeak \
  --Re2_BulkPeak /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/"${Re2[SLURM_ARRAY_TASK_ID]}"/macs2_temp/bulk_peaks_peaks.narrowPeak \
  --Re1_AllReads /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Re1[SLURM_ARRAY_TASK_ID]}"_Unique.bed \
  --Re2_AllReads /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Re2[SLURM_ARRAY_TASK_ID]}"_Unique.bed \
  --OutFileName "${OutFileNames[SLURM_ARRAY_TASK_ID]}" \
  --OutPath /scratch/sb14489/3.scATAC/2.Maize_ear/9.CheckQC/3.Replicates_Corr/
