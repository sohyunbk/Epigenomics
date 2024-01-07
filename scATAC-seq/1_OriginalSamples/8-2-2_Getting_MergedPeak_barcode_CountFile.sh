#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11.%j.err    # Standard error log
#SBATCH --array=0-27

ClusterN=(L1atFloralMeristem  Unknown1 FloralMeristem_SuppressedBract  PhloemPrecursor  Unknown2
G2_M ProcambialMeristem_ProtoXylem_MetaXylem  Unknown_Sclerenchyma IM-OC ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
Unknown_lowFRiP L1 SPM-base_SM-base  XylemParenchyma_PithParenchyma
L1atFloralMeristem  Unknown1 FloralMeristem_SuppressedBract  PhloemPrecursor  Unknown2
G2_M ProcambialMeristem_ProtoXylem_MetaXylem  Unknown_Sclerenchyma IM-OC ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
Unknown_lowFRiP L1 SPM-base_SM-base  XylemParenchyma_PithParenchyma)

ml Anaconda3/2023.09-0
source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Combined_Sorted_k12.bed_OnlyChr  \
  -b /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_MergedDifferentSizePeak/A619Bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_MergedPeak_Intergenic.bed  \
  -wa -wb \
  -sorted  | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl \
   - > /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_PeaksCount_byA619Barcode.txt

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Combined_Sorted_k12.bed_OnlyChr  \
 -b /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_MergedDifferentSizePeak/A619Bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_MergedPeak_Intergenic.bed  \
 -wa -wb \
 -sorted  | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl \
  - > /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_PeaksCount_byBif3Barcode.txt
