#!/bin/bash
#SBATCH --job-name=GetGeneMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/GetGeneMatrix.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/GetGeneMatrix.%j.err    # Standard error log
#SBATCH --array=0-1
ml Anaconda3/2023.09-0
conda activate r_env
SampleName=(A619 Bif3)
MetaFiles=(A619/Ref_AnnV4_metadata.txt Bif3/Bif3_AnnV4_metadata.txt)

Rscript /home/sb14489/Epigenomics/scATAC-seq/1_OriginalSamples/11-2-1_GetGeneXCTTable.sh \
--GA /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_"${SampleName[SLURM_ARRAY_TASK_ID]}"_includingZmCLE7.txt \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${MetaFiles[SLURM_ARRAY_TASK_ID]}" \
--OutputDir /scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC \
--SampleName "${SampleName[SLURM_ARRAY_TASK_ID]}"_AnnV4
