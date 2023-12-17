#!/bin/bash
#SBATCH --job-name=MarkerGene        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-5

ml Anaconda3/2022.10
source activate r_env

Sparse_Dir=(rel2_includingZmCLE7
rel2_includingZmCLE7
rel2_includingZmCLE7
Bif3_IncludingZmCLE7
A619_IncludingZmCLE7
A619_IncludingZmCLE7)
Metas=(rel2/rel2_Cluster2_Recluster_Sub_res1_knear100_Partmetadata.txt
rel2/rel2_Cluster3_Recluster_Sub_res1_knear100_Partmetadata.txt
rel2/rel2_Cluster4_Recluster_Sub_res1_knear100_Partmetadata.txt
Bif3/bif3_Cluster1_Recluster_Sub_res1_knear100_Partmetadata.txt
A619_PreviousScript/Cluster1_Recluster_Sub_res1_knear100_Partmetadata.txt
A619_PreviousScript/Cluster3_Recluster_Sub_res1_knear100_Partmetadata.txt)
OutDir=(rel2_ReCluster
rel2_ReCluster
rel2_ReCluster
Bif3_ReCluster
A619_ReCluster
A619_ReCluster)
OutFileNames=(rel2_SubCluster2_InsituMarkerGenes
rel2_SubCluster3_InsituMarkerGenes
rel2_SubCluster4_InsituMarkerGenes
Bif3_SubCluster1_InsituMarkerGenes
A619_SubCluster1_InsituMarkerGenes
A619_SubCluster3_InsituMarkerGenes)


Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/MarkerGenes_UMAPVisual_ForSubCluster.R \
 --imputed_sparse /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/"${Sparse_Dir[SLURM_ARRAY_TASK_ID]}"/opt_allgenes_impute.activity.rds \
 --meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/"${Metas[SLURM_ARRAY_TASK_ID]}" \
 --gene /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt \
 --OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/5.Subcluster_markergene/"${OutDir[SLURM_ARRAY_TASK_ID]}"/ \
 --prefix "${OutFileNames[SLURM_ARRAY_TASK_ID]}"


#  --gene /scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant.txt \
#  --OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/5.Subcluster_markergene/rel2_SubCluster/ \
#  --prefix rel2_SubCluster"${Cluesters[SLURM_ARRAY_TASK_ID]}"_DenovoA619Genes
