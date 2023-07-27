#!/bin/bash
#SBATCH --job-name=QC_JBrowser        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=2             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-2_QC_JBrowser.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-2_QC_JBrowser.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-3                   # Array range

List=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)

ml Anaconda3/2020.02
source activate r_env

DataDir="/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/"
BedDir="/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"

## 1) Filter_Bedfile_onlyfor high quality cells
#fgrep -f "$DataDir""${List[SLURM_ARRAY_TASK_ID]}"/"${List[SLURM_ARRAY_TASK_ID]}"_FilteredCellBarcode.txt \
#  "$BedDir""${List[SLURM_ARRAY_TASK_ID]}""_Unique.bed" > \
#  "$DataDir"/"${List[SLURM_ARRAY_TASK_ID]}"/Filtered_"${List[SLURM_ARRAY_TASK_ID]}".bed

## 2) Macs2
#module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4
cd "$DataDir"/"${List[SLURM_ARRAY_TASK_ID]}"/

#macs2 callpeak -t Filtered_"${List[SLURM_ARRAY_TASK_ID]}".bed -f BED --nomodel \
#                --keep-dup all --extsize 150 --shift -50 --qvalue .05 --bdg \
#                -n Filtered"${List[SLURM_ARRAY_TASK_ID]}"_macs2


### 3) bdg --> bw

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/Jbrowse

module load  SAMtools/1.10-iccifort-2019.5.281
module load BEDTools/2.29.2-GCC-8.3.0

sh /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.sh \
  -Step bdgTobw -bdgFile Filtered"${List[SLURM_ARRAY_TASK_ID]}"_macs2_treat_pileup.bdg \
 -Fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
  -OutputName Filtered"${List[SLURM_ARRAY_TASK_ID]}"_macs2_treat_pileup
