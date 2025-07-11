#!/bin/bash
#SBATCH --job-name=FixingBarcode        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcodes.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcodes.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-3                   # Array range

List=(1_A619 1_A619_2 3_bif3 3_bif3_2)
Path=/scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3/

cd "$Path"
mkdir -p "$Path"/4.Bam_FixingBarcode

module load Anaconda3/2022.10
source activate /home/sb14489/.conda/envs/ucsc
module load  SAMtools/1.10-GCC-8.3.0

#FixingBarcode
~/.conda/envs/ucsc/bin/python ../workflow_scripts/SocratesStart_QC/FixingBarcodeName.py \
 -BAM ./3.SortedBam/"${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam -exp_name "${List[SLURM_ARRAY_TASK_ID]}" | samtools view -@ 12 - > ./4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_BarcodeFixed.sam

 #FixingBarcode
 ~/.conda/envs/ucsc/bin/python ../workflow_scripts/SocratesStart_QC/MakeTn5bed.py \
 -sam ./4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_BarcodeFixed.sam -output_file ./4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_Unique.bed
