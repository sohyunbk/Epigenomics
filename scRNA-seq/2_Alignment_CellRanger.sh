#!/bin/bash
#SBATCH --job-name=Alignment_lowTime        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=28             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request
#SBATCH --time=80:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --array=0-3                   # Array range

List=(Sohyun-wt-1 Sohyun-wt-2 Sohyun-bif3-1 Sohyun-bif3-2)

Path=/scratch/sb14489/4.scRNAseq/2.snRNA-seq/
Reference=/scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNACellRanger_ZmCLE7/

module load CellRanger/7.0.0

cd "$Path"/2.Mapped

cellranger count --id="${List[SLURM_ARRAY_TASK_ID]}" \
                      --transcriptome=$Reference \
                      --fastqs=$Path/1.Raw_Data/ \
                      --sample="${List[SLURM_ARRAY_TASK_ID]}" \
                      --localcores=28 \
                      --localmem=300
