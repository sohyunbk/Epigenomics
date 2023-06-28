#!/bin/bash
#SBATCH --job-name=9-5        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9-5.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9-5.%j.err    # Standard error log
#SBATCH --array=0-1

SpareDdataName=(1_A619_Combined_Sorted_k12.bed 3_bif3_Combined_Sorted_k12.bed)
MetaData=(/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt)
OutputNames=(1_A619_Combined_AnnV3.bed 3_bif3_Combined_AnnV3.bed)

python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/7_PeakCalling_Browser/7-5_FilterSparseFileForFurtherStudy.py \
 -SparseBed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${SpareDdataName[SLURM_ARRAY_TASK_ID]}" \
 -MetaFile "${MetaData[SLURM_ARRAY_TASK_ID]}" \
 -OutputName /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${OutputNames[SLURM_ARRAY_TASK_ID]}"
