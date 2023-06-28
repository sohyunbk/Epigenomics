#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9_PeakCalling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9_PeakCalling.%j.err    # Standard error log
#SBATCH --array=0-3

Path=/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V2
Combined=(1_A619 2_rel2 3_bif3 4_relk1)
Sample=(A619 rel2 bif3 relk1)

ml Anaconda3/2020.02
source activate r_env

cd $Path

#~/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/9-2_call_scACRs.py -bed /scratch/sb14489/3.scATAC_flo/4.Bam_FixingBarcode/"${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted.bed -meta /scratch/sb14489/3.scATAC_flo/5.Socrates/4_CombineAll_AfterD/Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.Clustering_BasedOnRefAllCells.metadata.ClusterName.txt -col Ref_Cluster -base "${Sample[SLURM_ARRAY_TASK_ID]}" -outdir "${Sample[SLURM_ARRAY_TASK_ID]}"_PeakCalling -fai /scratch/sb14489/3.scATAC_flo/0.Reference/genome.fa.fai -bw TRUE

~/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/7_PeakCalling_Browser/7-2_call_scACRs.py \
 -bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted.bed \
 -meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/AnnotationV2."${Sample[SLURM_ARRAY_TASK_ID]}".metadata.txt \
 -col Ann_V2 -base "${Sample[SLURM_ARRAY_TASK_ID]}" -outdir "${Sample[SLURM_ARRAY_TASK_ID]}" \
 -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai -bw TRUE
