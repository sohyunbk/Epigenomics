#!/bin/bash
#SBATCH --job-name=Meme_motif        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Meme_motif.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Meme_motif.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0

BedFileName=(IM-OC.FDR0.01Bif3Higher IM-OC.FDR0.01A619Higher IM-OC.FDR0.05A619Higher IM-OC.FDR0.05Bif3Higher)

module load MEME/5.5.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.3.0

## Step1) Make null distribution
~/.conda/envs/r_env/bin/python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Meme/Generate_null_bedsample_forSTREAM.py \
--bed_file /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${BedFileName[SLURM_ARRAY_TASK_ID]}".bed \
--genome_file /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
--genome_index /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai \
--AllPeakForControl /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_MergedDifferentSizePeak/A619Bif3_IM-OC_MergedPeak_Intergenic.bed \
--output_name /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/"${BedFileName[SLURM_ARRAY_TASK_ID]}".ControlfromIntergenicAllSameCTPeaks.fa \
--Region Within
## Step2) STREME
bash /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Meme/Motif_Meme_FromACRBed.sh \
--infile_Bed /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${BedFileName[SLURM_ARRAY_TASK_ID]}".bed \
--Fa /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
--MemeMotifDB /scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt \
--ControlFA /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/"${BedFileName[SLURM_ARRAY_TASK_ID]}".ControlfromIntergenicAllSameCTPeaks.fa \
--OutfilePath /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/"${BedFileName[SLURM_ARRAY_TASK_ID]}".ControlfromIntergenicAllSameCTPeaks.XSTREME
