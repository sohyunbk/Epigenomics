#!/bin/bash
#SBATCH --job-name=LowMemory        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC/0.log/ChangeFile.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC/0.log/ChangeFile.%j.err    # Standard error log
#SBATCH --array=0-7                   # Array range

## This ExampleFailed As it has only two peaks
SampleName=(1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2)

cd /scratch/sb14489/3.scATAC/1.MaizeExample/1.Rawdata/OriginalData
List=`find -name "${SampleName[SLURM_ARRAY_TASK_ID]}""*" | sed 's|.*/||'`

for i in $List;
do
#zcat $Data_Path$i | head -n 400000 > $Output_Path/$i
split -a 4 -d -l $i $i__
done
