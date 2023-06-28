#!/bin/bash
#SBATCH --job-name=LowMemory        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=50gb                     # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/CuttingFiles.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/CuttingFiles.%j.err    # Standard error log
#SBATCH --array=0-7                   # Array range

## This ExampleFailed As it has only two peaks
SampleName=(1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2)
python /home/sb14489/1.scATAC-seq/1_scATAC-seq/1_MaizeExample/2_MakeExampleinBarcodeCount/CuttingUniqueBarcodeFile.py "${SampleName[SLURM_ARRAY_TASK_ID]}"
