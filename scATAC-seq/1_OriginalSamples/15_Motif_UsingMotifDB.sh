#!/bin/bash
#SBATCH --job-name=Meme_motif        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Meme_motif.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Meme_motif.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

module load MEME/5.5.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.3.0

#Infile_Bed="/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.FDR0.05.Bed"
#Fafile="/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa"
#MemeMotifDB="/scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt"
#OutfilePathName="/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_JASPARMotif"

bash /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/FindMotif_FromPeak.sh \
--Infile_Bed /scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.A619Higher.Bed \
--Fafile /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
--MemeMotifDB /scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt \
--OutfilePathName /scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_A619Higher_JASPARMotif \

bash /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/FindMotif_FromPeak.sh \
--Infile_Bed /scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC.Bif3Higher.Bed \
--Fafile /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa \
--MemeMotifDB /scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt \
--OutfilePathName /scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_Bif3Higher_JASPARMotif \
