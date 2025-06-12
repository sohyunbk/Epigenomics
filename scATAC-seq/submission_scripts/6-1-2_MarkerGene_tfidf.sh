#!/bin/bash
#SBATCH --job-name=MarkerGene_rel2        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=30gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-7

source activate r_env

WDNames=(rel2_includingZmCLE7
rel2_TopDenovoMarkersA619
Bif3_IncludingZmCLE7
Bif3_TopDenovoMarkersA619
A619_IncludingZmCLE7
A619_TopDenovoMarkersA619
A619_WOX_ARF
Bif3_WOX_ARF)
FileNames=(rel2_includingZmCLE7
rel2_TopDenovoMarkersA619
Bif3_IncludingZmCLE7
Bif3_TopDenovoMarkersA619
A619_IncludingZmCLE7
A619_TopDenovoMarkersA619
A619_WOX_ARF
Bif3_WOX_ARF)
Metas=(6.Annotation/0.AnnotatedMeta/rel2/rel2_AnnV4.txt
6.Annotation/0.AnnotatedMeta/rel2/rel2_AnnV4.txt
6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV3_metadata.txt
6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV3_metadata.txt
6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt
6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt
6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt
6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV3_metadata.txt)
GAs=(GA_rel2_includingZmCLE7.txt
GA_rel2_includingZmCLE7.txt
GA_Bif3_Re.txt
GA_Bif3_Re.txt
GA_A619_Re.txt
GA_A619_Re.txt
GA_A619_Re.txt
GA_Bif3_Re.txt)
PCAs=(Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt
Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt
Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt
Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
)
Markers=(230426_EarMarker.txt
231113_Top5DenovoGenesinA619_NoRedundant.txt
230426_EarMarker.txt
231113_Top5DenovoGenesinA619_NoRedundant.txt
230426_EarMarker.txt
231113_Top5DenovoGenesinA619_NoRedundant.txt
ARF_WOX.txt
ARF_WOX.txt)

Rscript ../workflow_scripts/Annotation_Cluster/MarkerGenes_UMAPVisual.R \
--WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/"${WDNames[SLURM_ARRAY_TASK_ID]}" \
--Name "${FileNames[SLURM_ARRAY_TASK_ID]}" \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/"${Metas[SLURM_ARRAY_TASK_ID]}" \
--geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${GAs[SLURM_ARRAY_TASK_ID]}" \
--pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/"${PCAs[SLURM_ARRAY_TASK_ID]}" \
--markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/"${Markers[SLURM_ARRAY_TASK_ID]}"
