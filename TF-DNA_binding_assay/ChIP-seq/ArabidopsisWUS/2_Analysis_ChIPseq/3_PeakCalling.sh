#!/bin/bash
#SBATCH --job-name=macs2        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/macs2.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/macs2.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

#bowtie2-build --threads 20 /scratch/sb14489/0.Reference/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /scratch/sb14489/0.Reference/TAIR10/TAIR10
mkdir /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling

macs2 callpeak \
 -t /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/2.Mapped/SRR8192660_unique_bowtie2_algn.bam \
	-c /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/2.Mapped/SRR8192661_unique_bowtie2_algn.bam \
 	-f BAM -g 1.3e+8 --nomodel   --extsize 100 -q  0.01 \
	-n WUS_GS_NonModel \
	--outdir /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling/

#134,634,692
#1000,000,000
#1.0e+9

#10833
