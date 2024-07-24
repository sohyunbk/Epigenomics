#!/bin/bash
#SBATCH --job-name=DownloadSRA       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=80gb                     # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC_flo/100.log/DownloadSRA.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC_flo/100.log/DownloadSRA.%j.err    # Standard error log
#SBATCH --array=0-2                   # Array range

module load SRA-Toolkit/2.11.1-centos_linux64

Files=(SRR12276825 SRR12276805 SRR12276824)
NewName=(B73re1_S1_L001 B73re2_S2_L002 B73re3_S3_L003)
NewName2=(B73re1 B73re2 B73re3)

cd /scratch/sb14489/4.scRNAseq/1.Raw_Data

## 1) Download big sra files
#prefetch --max-size 40g "${Files[SLURM_ARRAY_TASK_ID]}"

#SRR12276825 Re1 -  PRJNA646989  - 39.4G bases, 24.4Gb downloads
#SRR12276805 Re2 - PRJNA646996  -  28.6G bases, 13.5Gb downloads
#SRR12276824 Re3  - PRJNA647001  - 42G bases, 24.5Gb downloads

## 2) Convert sra to fastq
#scRNA-seq libraries were sequenced by Illumina short reads with 400M paired end reads per library (read1 = 28bp, read2 = 56bp
#fasterq-dump -S --include-technical "${Files[SLURM_ARRAY_TASK_ID]}"/"${Files[SLURM_ARRAY_TASK_ID]}".sra

## 3) Change file names
#SRR12276805_1.fastq, SRR12276805_2.fastq, SRR12276805_3.fastq
#1_A619_S9_L002_I1_001.fastq.gz  1_A619_S9_L002_R1_001.fastq.gz  1_A619_S9_L002_R2_001.fastq.gz  1_A619_S9_L002_R3_001.fastq.gz
#[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# _1 : 8bp --> I1
# _2: 28bp --> R1 (16bp barcode + 10bp UMI)
# _3: 56bp --> R2
#mv "${Files[SLURM_ARRAY_TASK_ID]}"_1.fastq "${NewName[SLURM_ARRAY_TASK_ID]}"_I1_001.fastq
#mv "${Files[SLURM_ARRAY_TASK_ID]}"_2.fastq "${NewName[SLURM_ARRAY_TASK_ID]}"_R1_001.fastq
#mv "${Files[SLURM_ARRAY_TASK_ID]}"_3.fastq "${NewName[SLURM_ARRAY_TASK_ID]}"_R2_001.fastq

mkdir  "${NewName2[SLURM_ARRAY_TASK_ID]}"
mv "${NewName2[SLURM_ARRAY_TASK_ID]}"_* "${NewName2[SLURM_ARRAY_TASK_ID]}"
