#!/bin/bash
#SBATCH --job-name=SubCluster        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/SubCluster.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/SubCluster.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-3

ml Anaconda3/2020.02
source activate r_env

Cluesters=(1 3 4)

module load CellRanger-ATAC/2.0.0

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/Subclustering.R \
--SampleName rel2 \
--MetaFile /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt \
--ObjAfterHarmony /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.rds \
--AnnSlot LouvainClusters \
--TargetClusterName "${Cluesters[SLURM_ARRAY_TASK_ID]}" \
--OutputDir /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/rel2
