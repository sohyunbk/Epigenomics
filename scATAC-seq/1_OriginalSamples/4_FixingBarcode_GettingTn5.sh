#!/bin/bash
#SBATCH --job-name=FixingBarcode_GetTn5        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcode_GetTn5.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcode_GetTn5.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0

List=(2_rel2_2)

## Idk why but 2_rel2_2_Unique.bed was missing

module load Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc

module load  SAMtools/1.10-iccifort-2019.5.281

#samtools index -@ 20 /scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/"${SampleNameList[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam

#FixingBarcode
#python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/FixingBarcodeName.py \
# -BAM /scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/"${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam \
# -exp_name "${List[SLURM_ARRAY_TASK_ID]}" | samtools view -@ 12 - > /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_BarcodeFixed.sam

 #FixingBarcode
 python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MakeTn5bed_fromBam.py \
 -sam /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_BarcodeFixed.sam \
 -output_file /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_Unique.bed
