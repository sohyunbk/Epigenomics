#!/bin/bash
#SBATCH --job-name=RemoveMultiMap_Deduplication_Browserr        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=400gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/RemoveMultiMap_Deduplication_Browserr.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/RemoveMultiMap_Deduplication_Browserr.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

OGSampleNameList=(Sohyun_A619-1 Sohyun_A619-2 Sohyun_BIF3-1 Sohyun_BIF3-2)
NewSampleNameList=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)

module load Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/r_env
module load picard/2.16.0-Java-1.8.0_144
module load  SAMtools/1.10-iccifort-2019.5.281
module load BEDTools/2.29.2-GCC-8.3.0

## RemovemultiMap & deduplication
sh ../workflow_scripts/RemoveMultiMap_Deduplication.sh --path /scratch/sb14489/3.scATAC/2.Maize_ear/ \
--MappedDir 2.Mapped_CellRangerv2  --OGSampleName "${OGSampleNameList[SLURM_ARRAY_TASK_ID]}" \
 --NewSampleName_forBam "${NewSampleNameList[SLURM_ARRAY_TASK_ID]}"

## GenomeBrowser
ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc
module load SAMtools/1.10-iccifort-2019.5.281
module load BEDTools/2.29.2-GCC-8.3.0

bash /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.sh -Step BamTobw  \
 -Fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
  -bam /scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/"${NewSampleNameList[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam
