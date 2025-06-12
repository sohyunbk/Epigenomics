#!/bin/bash
#SBATCH --job-name=QC        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-1_QC.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-1_QC.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-3                   # Array range

List=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)

ml Anaconda3/2022.10
source activate r_env

Rscript ../workflow_scripts/MakeSocratesObject_QCFilterCells.R  \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55_PeakClustering/ \
 --BinSize peak --Name "${List[SLURM_ARRAY_TASK_ID]}" \
 --bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_Unique.bed \
 --ann /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf \
 --chr /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
 --MinTn5 1000 --TSS 0.35 --FRiP 0.55 --Org 0.1 --Step All
