#!/bin/bash
#SBATCH --job-name=Alignment_lowTime        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=100             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/BuildRef.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/BuildRef.%j.err    # Standard error log

Path=/scratch/sb14489/0.Reference/Maize_B73/
Reference=Maize_B73_V5_withMtPt_scRNACellRanger_ZmCLE7

module load CellRanger/7.0.0

cd "$Path"

cellranger mkref \
  --genome="$Reference" \
  --fasta=Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
  --genes=Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf_AddZmCLE7.gtf --nthreads=100 --memgb=100
