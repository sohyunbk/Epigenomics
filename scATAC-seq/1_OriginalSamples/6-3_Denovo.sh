#!/bin/bash
#SBATCH --job-name=Denovo        # Job name
#SBATCH --partition=sbatch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Denovo.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Denovo.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail


ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/De_novo_marker_byCluster.R \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt \
--GeneBA /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619.txt \
--marker /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt \
--bed /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed \
--Ann_ColumnName LouvainClusters \
--Species maize --OutputBaseName Ref_LouvainClusters \
--OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV3
