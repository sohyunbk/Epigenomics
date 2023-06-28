#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11_Intersect.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11_Intersect.%j.err    # Standard error log
#SBATCH --array=0-12

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

#conda source /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_For_dACR

bedtools intersect -wo -a "${ClusterN[SLURM_ARRAY_TASK_ID]}"_A619_bif3_LongerPeak_Sorted.bed \
  -b /scratch/sb14489/0.Reference/Maize_B73/Zm.final_blaclist.final.txt   > \
  "${ClusterN[SLURM_ARRAY_TASK_ID]}"_Overlap_ComA619Bif3_BlackList.bed

bedtools intersect -wo -a "${ClusterN[SLURM_ARRAY_TASK_ID]}"_A619_bif3_LongerPeak_Sorted.bed \
  -b /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed   > \
  "${ClusterN[SLURM_ARRAY_TASK_ID]}"_Overlap_ComA619Bif3_Genes.bed
