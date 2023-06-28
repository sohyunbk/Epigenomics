#!/bin/bash
#SBATCH --job-name=LowMemory        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC/0.log/ChangeFile.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC/0.log/ChangeFile.%j.err    # Standard error log
#SBATCH --array=0-7                   # Array range

## This ExampleFailed As it has only two peaks 
SampleName=(1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2)
SampleName_Change=(ExA619Re1 ExRel2Re1  ExBif3Re1 ExRelk1Re1 ExA619Re2 ExRel2Re2 ExBif3Re2 ExRelk1Re2)

Data_Path=/scratch/sb14489/3.scATAC_flo_Past/1.Rawdata/"${SampleName[SLURM_ARRAY_TASK_ID]}"/
RawDataFormat=".fastq.gz"
Output_Path=/scratch/sb14489/3.scATAC/1.MaizeExample/1.Rawdata
List=`find "$Data_Path" -name "*""$RawDataFormat" | sed 's|.*/||'`

for i in $List;
do
#echo $item
#<<"END"
ChangedName_v1=${i/$SampleName/"${SampleName_Change[SLURM_ARRAY_TASK_ID]}"}
ChangedName=${ChangedName_v1/".gz"/""}

#zcat $Data_Path$i | head -n 400000 > $Output_Path/$i
echo "$Output_Path/"${SampleName_Change[SLURM_ARRAY_TASK_ID]}""
echo "zcat $Data_Path$i | head -n 400000 > $Output_Path/"${SampleName_Change[SLURM_ARRAY_TASK_ID]}"/$ChangedName"
mkdir $Output_Path/"${SampleName_Change[SLURM_ARRAY_TASK_ID]}"
zcat $Data_Path$i | head -n 400000 > $Output_Path/"${SampleName_Change[SLURM_ARRAY_TASK_ID]}"/$ChangedName

done
