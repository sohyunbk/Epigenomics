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

ml Anaconda3/2022.10
source activate r_env

Cluesters=(1 3 4)

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/MarkerGenes_UMAPVisual_ForSubCluster.R \
 --imputed_sparse /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/rel2_includingZmCLE7/opt_allgenes_impute.activity.rds \
 --meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/rel2/rel2_Cluster"${Cluesters[SLURM_ARRAY_TASK_ID]}"_Recluster_Sub_res1_knear100_Partmetadata.txt \
 --gene /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt \
 --OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/5.Subcluster_markergene/rel2_SubCluster/ \
 --prefix rel2_SubCluster"${Cluesters[SLURM_ARRAY_TASK_ID]}"_InsituMarkerGenes


 Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/MarkerGenes_UMAPVisual_ForSubCluster.R \
  --imputed_sparse /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/rel2_includingZmCLE7/opt_allgenes_impute.activity.rds \
  --meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/rel2/rel2_Cluster"${Cluesters[SLURM_ARRAY_TASK_ID]}"_Recluster_Sub_res1_knear100_Partmetadata.txt \
  --gene /scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant.txt \
  --OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/5.Subcluster_markergene/rel2_SubCluster/ \
  --prefix rel2_SubCluster"${Cluesters[SLURM_ARRAY_TASK_ID]}"_DenovoA619Genes
