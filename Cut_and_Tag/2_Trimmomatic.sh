#!/bin/bash
#SBATCH --job-name=Trimmomatic        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=30gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Trimmomatic.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Trimmomatic.%j.err    # Standard error log
#SBATCH --array=0-7

module load Trimmomatic/0.39-Java-13

cd /scratch/sb14489/7.DAPorChIP/CUTandTAG/

SampleNames=(SRR21185825 SRR21228818 SRR21228819 SRR21228820 SRR21228821 SRR21228813 SRR21228814 SRR21228816 SRR21228817)
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 20 \
 ./1.RawData/"${SampleNames[SLURM_ARRAY_TASK_ID]}"_1.fastq ./1.RawData/"${SampleNames[SLURM_ARRAY_TASK_ID]}"_2.fastq \
 ./2.Trimmomatic/"${SampleNames[SLURM_ARRAY_TASK_ID]}"_forward_paired.fastq ./2.Trimmomatic/"${SampleNames[SLURM_ARRAY_TASK_ID]}"_forward_unpaired.fastq \
 ./2.Trimmomatic/"${SampleNames[SLURM_ARRAY_TASK_ID]}"_reverse_paired.fastq ./2.Trimmomatic/"${SampleNames[SLURM_ARRAY_TASK_ID]}"_reverse_unpaired.fastq \
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
