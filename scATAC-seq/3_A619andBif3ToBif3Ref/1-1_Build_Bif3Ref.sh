#!/bin/bash
#SBATCH --job-name=BuildReference       # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=80gb                     # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/BuildRef.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/BuildRef.%j.err    # Standard error log

Path=/scratch/sb14489/0.Reference/Maize_Ki3

cd $Path
module load CellRanger-ATAC/2.0.0
cellranger-atac mkref --config=CellRangerv2-ATAC_Ki3_Bif3.config ## In protocl they are saying it should be gtf but gff seems fine in v1 but in v2 it should be gtf
