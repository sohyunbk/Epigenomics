#!/bin/bash
#SBATCH --job-name=3_QC        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=6:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-3_CleanFeatureDetectDoublet.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-3_CleanFeatureDetectDoublet.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-3                   # Array range

List=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/CleanFeatures_DetectDoublets.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/ \
 --Name "${List[SLURM_ARRAY_TASK_ID]}" --PreFix Tn5Cut1000_Binsize500 --MinT 0.007 --MaxT 0.005 --nPC 100
# --Name "${List[SLURM_ARRAY_TASK_ID]}" --PreFix Tn5Cut1000_Binsize500 --MinT 0.001 --MaxT 0.05 --nPC 100
# --Name "${List[SLURM_ARRAY_TASK_ID]}" --PreFix Tn5Cut1000_Binsize500 --MinT 0.01 --MaxT 0.05 --nPC 100

#Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/CleanFeatures_DetectDoublets.R \
# --WD /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/ \
# --Name "${List[SLURM_ARRAY_TASK_ID]}" --PreFix Tn5Cut1000_Binsize500 --MinT 0.01 --MaxT 0.05 --nPC 100
