#!/bin/bash
#SBATCH --job-name=MarkerGene_rel2        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
--WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/rel2_includingZmCLE7 \
--Name rel2_ZmCLE7 \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt \
--geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_rel2_includingZmCLE7.txt \
--pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
--markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
--WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/rel2_TopDenovoMarkersA619 \
--Name rel2_TopDenovoMarkersA619 \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt \
--geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_rel2_includingZmCLE7.txt \
--pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
--markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant.txt

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
--WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_IncludingZmCLE7_ReTry \
--Name Bif3_ZmCLE7 \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.metadata.txt \
--geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_Bif3_Re.txt \
--pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
--markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
--WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_TopDenovoMarkersA619 \
--Name Bif3_TopDenovoMarkersA619 \
--meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.metadata.txt \
--geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_Bif3_Re.txt \
--pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
--markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant.txt


Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_TopDenovoMarkersA619 \
 --Name A619_TopDenovoMarkersA619 \
 --meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV4_metadata.txt \
 --geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re.txt \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt \
 --markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant.txt
