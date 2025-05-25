#!/bin/bash
#SBATCH --job-name=Fimo        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=7:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Fimo.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Fimo.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-12
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

Celltype=(FloralMeristem_SuppressedBract G2_M IM-OC L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base Unknown1 Unknown2 Unknown_Sclerenchyma Unknown_lowFRiP XylemParenchyma_PithParenchyma)

module load MEME/5.5.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.3.0

MemeMotifDB=/scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/CentralZone_TAAT.txt
Infile_bed=/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${Celltype[SLURM_ARRAY_TASK_ID]}".FDR0.05Bif3Higher.bed
OutfilePathName=/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/3.fimo/Revision_"${Celltype[SLURM_ARRAY_TASK_ID]}".FDR0.05Bif3Higher
Infile_fasta="${Infile_bed%.bed}.fasta"

bedtools getfasta -fi /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa -bed "$Infile_bed" -fo "$Infile_fasta"
fimo --o "$OutfilePathName" "$MemeMotifDB" "$Infile_fasta"
