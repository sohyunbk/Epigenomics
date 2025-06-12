#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/dACR_dotFigure.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/dACR_dotFigure.%j.err    # Standard error log

ml Anaconda3/2023.09-0
source activate r_env

Rscript ../workflow_scripts/Viualization/dACR_number_byCellType_DotPlot_withEdgeRSummaryFile.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/ \
 --FDRCutOff 0.01 \
 --OutFilename dACRNumber_DotPlot_withoutTotal_FDR0.01_RemoveUnkown_FixingAxis.pdf \
 --CellOrderFile /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3_Reverse.txt

 Rscript ../workflow_scripts/Viualization/dACR_number_byCellType_DotPlot_withEdgeRSummaryFile.R \
  --WD /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/ \
  --FDRCutOff 0.05 \
  --OutFilename dACRNumber_DotPlot_withoutTotal_FDR0.05_RemoveUnkown_FixingAxis.pdf \
  --CellOrderFile /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3_Reverse.txt
