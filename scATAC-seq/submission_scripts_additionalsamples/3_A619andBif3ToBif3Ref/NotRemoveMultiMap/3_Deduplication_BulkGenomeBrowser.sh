#!/bin/bash
#SBATCH --job-name=RemoveMultiMap_Deduplication_Browserr        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=400gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=50:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/RemoveMultiMap_Deduplication_Browserr.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/RemoveMultiMap_Deduplication_Browserr.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

OGSampleNameList=(1_A619  1_A619_2  3_bif3  3_bif3_2)
NewSampleNameList=(1_A619  1_A619_2  3_bif3  3_bif3_2)

module load Anaconda3/2022.10
source activate /home/sb14489/.conda/envs/r_env
module load picard/2.25.1-Java-11
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.29.2-GCC-8.3.0

## RemovemultiMap & deduplication
sh ../workflow_scripts/Mapping_RefiningBam/RemoveMultiMap_Deduplication_PicardVersionUpdate.sh \
 --path /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3_NotRemoveMultiMap/ --RemoveDup No \
--MappedDir 2.Mapped  --OGSampleName "${OGSampleNameList[SLURM_ARRAY_TASK_ID]}" \
 --NewSampleName_forBam "${NewSampleNameList[SLURM_ARRAY_TASK_ID]}"

## GenomeBrowser
module load Anaconda3/2022.10
source activate /home/sb14489/.conda/envs/ucsc
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.29.2-GCC-8.3.0

bash /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.sh -Step BamTobw  \
 -Fai /scratch/sb14489/0.Reference/Maize_Ki3/Zm-Ki3-REFERENCE-NAM-1.0_OnlyChr_Bif3.fa.fai \
  -bam /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3_NotRemoveMultiMap//3.SortedBam/"${NewSampleNameList[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam
