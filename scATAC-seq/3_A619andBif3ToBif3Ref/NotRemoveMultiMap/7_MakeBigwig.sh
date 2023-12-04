#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Bigwig.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Bigwig.%j.err    # Standard error log
#SBATCH --array=0-5
OutFileList=(1_A619_Re2 1_A619_Re1 3_bif3_Re2 3_bif3_Re1 A619_Re1andRe2 Bif3_Re1andRe2)

source activate /home/sb14489/.conda/envs/ucsc

# Directory path
directory_path=/scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3_NotRemoveMultiMap/5.Jbrowse_MACS2/"${OutFileList[SLURM_ARRAY_TASK_ID]}"
# Loop through files ending with ".macs_treat_pileup_CPM.bdg"
for file in "$directory_path"/*.macs_treat_pileup_CPM.bdg; do
    # Check if file exists to avoid processing an empty glob
    if [[ -f "$file" ]]; then
      echo "Processing file: $file"
      bash /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.sh \
      -Step bdgTobw -bdgFile $file \
      -Fai /scratch/sb14489/0.Reference/Maize_Ki3/Zm-Ki3-REFERENCE-NAM-1.0_OnlyChr_Bif3.fa.fai \
      -OutputName "${file%.bdg}.bw"
    fi
done
