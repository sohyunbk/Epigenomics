#!/bin/bash
#SBATCH --job-name=bowtie2        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/bowtie2.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/bowtie2.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-7

module load Bowtie2/2.5.2-GCC-11.3.0
Sample=(SRR21185825 SRR21228818 SRR21228819 SRR21228820 SRR21228821 SRR21228813 SRR21228814 SRR21228816 SRR21228817)

#bowtie2-build --threads 20 /scratch/sb14489/0.Reference/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /scratch/sb14489/0.Reference/TAIR10/TAIR10
cd /scratch/sb14489/7.DAPorChIP/CUTandTAG/


bowtie2 -x /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf \
-1 ./2.Trimmomatic/"${Sample[SLURM_ARRAY_TASK_ID]}"_forward_paired.fastq \
-2 ./2.Trimmomatic/"${Sample[SLURM_ARRAY_TASK_ID]}"_reverse_paired.fastq \
-b ./3.Bowtie2/"${Sample[SLURM_ARRAY_TASK_ID]}"_bowtie2_algn.bam -p 20 \
2> ./3.Bowtie2/"${Sample[SLURM_ARRAY_TASK_ID]}".log
