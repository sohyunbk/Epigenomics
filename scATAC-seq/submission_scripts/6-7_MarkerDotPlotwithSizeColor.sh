#!/bin/bash
#SBATCH --job-name=6-7_Dot        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6-7_Dot.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6-7_Dot.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-2

source activate r_env
GAFiles=(
GA_A619_Re.txt
GA_Bif3_Re.txt
GA_rel2_includingZmCLE7.txt
GA_A619_Re.txt
GA_Bif3_Re.txt
GA_rel2_includingZmCLE7.txt)
MetaFiles=(A619/Ref_AnnV4_metadata.txt
Bif3/Bif3_AnnV4_metadata.txt
rel2/rel2_AnnV4.txt
A619/Ref_AnnV4_metadata.txt
Bif3/Bif3_AnnV4_metadata.txt
rel2/rel2_AnnV4.txt)
MarkerGenes=(
SeletedMarkergeneForDotPlot_RemoveUnknown.txt
SeletedMarkergeneForDotPlot_RemoveUnknown.txt
SelectedMarkerGeneForDotPlot.txt
231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt
231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt
231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt)
OutPrefixs=(
A619_Annv4_InsituMarker_RemoveUnknown
Bif3_Annv4_InsituMarker_RemoveUnknown
Rel2_Annv4_InsituMarker
A619_Annv4_DenovoGenes
Bif3_Annv4_DenovoGenes
Rel2_Annv4_DenovoGenes)
CellOrderFiles=(Ann_v4_CellType_order_forA619Bif3_RemoveUnknown.txt
Ann_v4_CellType_order_forA619Bif3_RemoveUnknown.txt
Ann_v4_CellType_order_forrel2.txt
Ann_v4_CellType_order_forA619Bif3.txt
Ann_v4_CellType_order_forA619Bif3.txt
Ann_v4_CellType_order_forrel2.txt)
MarkerGeneOrderFile=(SelectedMarkerGeneForDotPot_GeneNameOrder_RemoveUnknown.txt
SelectedMarkerGeneForDotPot_GeneNameOrder_RemoveUnknown.txt
SelectedMarkerGeneForDotPot_GeneNameOrder.txt
231113_Top5DenovoGenesinA619_GeneNameOrder.txt
231113_Top5DenovoGenesinA619_GeneNameOrder.txt
231113_Top5DenovoGenesinA619_GeneNameOrder.txt)

Rscript ../workflow_scripts/Annotation_Cluster/DotPlot_MarkerGenes_withColorandDotSize.R \
--GA /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${GAFiles[SLURM_ARRAY_TASK_ID]}" \
--Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${MetaFiles[SLURM_ARRAY_TASK_ID]}" \
--MarkerGene /scratch/sb14489/3.scATAC/0.Data/MarkerGene/"${MarkerGenes[SLURM_ARRAY_TASK_ID]}" \
--OutPathandPrefix /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/7.DotPlot/"${OutPrefixs[SLURM_ARRAY_TASK_ID]}" \
--AnnSlot Ann_v4 \
--CellOrdertxt /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${CellOrderFiles[SLURM_ARRAY_TASK_ID]}" \
--MarkerOrdertxt /scratch/sb14489/3.scATAC/0.Data/MarkerGene/"${MarkerGeneOrderFile[SLURM_ARRAY_TASK_ID]}"
