#!/bin/bash
#SBATCH --job-name=CombineBedfiles        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=500gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9-1_CombineBeds_%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9-1_CombineBeds_%j.err    # Standard error log
#SBATCH --array=0-3

Re1=(1_A619 2_rel2 3_bif3 4_relk1)
Re2=(1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2) #Add Later #SBATCH --array=0-3
Combined=(1_A619 2_rel2 3_bif3 4_relk1)


cd /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode

module load BEDTools/2.30.0-GCC-12.2.0

cat "${Re1[SLURM_ARRAY_TASK_ID]}"_Unique.bed "${Re2[SLURM_ARRAY_TASK_ID]}"_Unique.bed  | grep "^chr"  > "${Combined[SLURM_ARRAY_TASK_ID]}"_Combined.bed
sort -k1,1 -k2,2n "${Combined[SLURM_ARRAY_TASK_ID]}"_Combined.bed > "${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted_k12.bed

## It's past sorting method
#bedtools sort -i "${Combined[SLURM_ARRAY_TASK_ID]}"_Combined.bed > "${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted.bed
## Bedtools sort requires lots of memories ........
