#!/bin/bash
#SBATCH --job-name=bowtie2        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=40:00:00               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/bowtie2.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/bowtie2.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-1                   # Array range

module load SAMtools/1.16.1-GCC-11.3.0

#bowtie2-build --threads 20 /scratch/sb14489/0.Reference/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /scratch/sb14489/0.Reference/TAIR10/TAIR10
Sample=(SRR8192660 SRR8192661)

grep -E "@|NM:" /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/2.Mapped/"${Sample[SLURM_ARRAY_TASK_ID]}"_bowtie2_algn.sam | \
grep -v "XS:" > /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/2.Mapped//"${Sample[SLURM_ARRAY_TASK_ID]}"_unique_bowtie2_algn.sam

samtools sort /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/2.Mapped//"${Sample[SLURM_ARRAY_TASK_ID]}"_unique_bowtie2_algn.sam \
 > /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/2.Mapped//"${Sample[SLURM_ARRAY_TASK_ID]}"_unique_bowtie2_algn.bam -@ 14
