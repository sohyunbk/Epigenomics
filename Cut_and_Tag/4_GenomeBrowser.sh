#!/bin/bash
#SBATCH --job-name=Jbrowse        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=2            # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=3:00:00               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/Jbrowse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Jbrowse.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-7

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/Jbrowse

module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-12.2.0

Sample=(SRR21185825 SRR21228818 SRR21228819 SRR21228820 SRR21228821 SRR21228813 SRR21228814 SRR21228816 SRR21228817)

#bowtie2-build --threads 20 /scratch/sb14489/0.Reference/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /scratch/sb14489/0.Reference/TAIR10/TAIR10
cd /scratch/sb14489/7.DAPorChIP/CUTandTAG/3.Bowtie2

bash /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.sh -Step BamTobw  \
 -Fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
  -bam "${Sample[SLURM_ARRAY_TASK_ID]}"_bowtie2_algn.bam
