#!/bin/bash
#SBATCH --job-name=Meme_motif        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Meme_motif.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Meme_motif.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

module load MEME/5.5.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.3.0

~/.conda/envs/r_env/bin/python ../workflow_scripts/Meme/Generate_null_bedsample_forSTREAM.py \
--bed_file /scratch/sb14489/3.scATAC/2.Maize_ear/14.CellTypeSpecificACRs/CellTypeACRs.A619_L1Specific_Intergenic.bed \
--genome_file /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
--genome_index /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai \
--AllPeakForControl /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic.bed \
--output_name /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic_ControlforL1Secific.fa \
--Region Within

bash ../workflow_scripts/Motif_Meme_FromACRBed.sh \
--infile_Bed /scratch/sb14489/3.scATAC/2.Maize_ear/14.CellTypeSpecificACRs/CellTypeACRs.A619_L1Specific_Intergenic.bed \
--Fa /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
--MemeMotifDB /scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt \
--ControlFA /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic_ControlforL1Secific.fa \
--OutfilePath /scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/CellTypeACRs.A619_L1Specific
