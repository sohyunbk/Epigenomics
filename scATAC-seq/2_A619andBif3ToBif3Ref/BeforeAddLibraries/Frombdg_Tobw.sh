#!/bin/bash
#SBATCH --job-name=ChangeFileForJbrowse        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=7:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ChangeFileForJbrowse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ChangeFileForJbrowse.%j.err    # Standard error log

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/Jbrowse

FileList=`find "/scratch/sb14489/3.scATAC/4.Bif3Ref/5.Jbrowse_MACS2" -name "*_pileup_CPM.bdg" | sed 's|.*/||'`

#for i in $FileList; do echo $i
#for i in $FileList; do echo "${i/_treat_pileup_CPM.bdg/}"

for i in $FileList;
do
python /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.py \
      -Step bdgTobw -bdgFile $i  -Fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai -OutputName "${i/_treat_pileup_CPM.bdg/}"
done
