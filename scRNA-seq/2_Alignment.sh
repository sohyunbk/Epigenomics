#!/bin/bash
#SBATCH --job-name=Alignment_lowTime        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request
#SBATCH --time=18:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --array=0-2                   # Array range

List=(B73re1 B73re2 B73re3)

Path=/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/
Reference=/scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNACellRanger/

module load CellRanger/7.0.0

cd "$Path"/2.Mapped

cellranger count --id="${List[SLURM_ARRAY_TASK_ID]}" \
                      --transcriptome=$Reference \
                      --fastqs=$Path/1.Raw_Data/ \
                      --sample="${List[SLURM_ARRAY_TASK_ID]}" \
                      --localcores=12 \
                      --localmem=300
