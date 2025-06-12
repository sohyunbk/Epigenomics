#!/bin/bash
#SBATCH --job-name=4_UMAP        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-4_UMAP_Comb2Re.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-4_UMAP_Comb2Re.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-1                   # Array range

List=(A619 bif3)

List1=(A619_Re3 bif3_Re3)
List2=(A619_Re4 bif3_Re4)

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/UMAP_2Replicates_UsingOtherWindow.R \
 --PreFix_name Tn5Cut1000_Binsize500_MinT0.005_MaxT0.01_PC100 \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS3_FRiP4// \
 --SampleS "${List[SLURM_ARRAY_TASK_ID]}"  \
  --Re1  "${List1[SLURM_ARRAY_TASK_ID]}" --Re2 "${List2[SLURM_ARRAY_TASK_ID]}"

  #Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/UMAP_2Replicates_UsingOtherWindow.R \
  # --PreFix_name Tn5Cut1000_Binsize500_MinT0.007_MaxT0.005_PC100 \
  # --WD /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55// \
  # --SampleS "${List[SLURM_ARRAY_TASK_ID]}"  \
  #  --Re1  "${List1[SLURM_ARRAY_TASK_ID]}" --Re2 "${List2[SLURM_ARRAY_TASK_ID]}"
