#!/bin/bash
#SBATCH --job-name=Barcode        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=500gb                     # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC/0.log/CuttingFiles.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC/0.log/CuttingFiles.%j.err    # Standard error log
#SBATCH --array=0-7                   # Array range

## This ExampleFailed As it has only two peaks
SampleName=(1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2)
python /home/sb14489/1.scATAC-seq/1_scATAC-seq/1_MaizeExample/1_BuildReference_ModifyFastq/PickUpCells_byBarcode_v2_faster.py "${SampleName[SLURM_ARRAY_TASK_ID]}"
