#!/bin/bash
#SBATCH --job-name=MarkerGene        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-2

source activate r_env

MetaFiles=(1 3 4)
Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/DotPlot_MarkerGenes_withColorandDotSize.R \
--GA /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re.txt \
--Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt \
--MarkerGene /scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt \
--OutPathandPrefix /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/7.DotPlot/A619_Annv4_DenovoGenes \
--AnnSlot Ann_v4 \
--CellOrdertxt /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt
