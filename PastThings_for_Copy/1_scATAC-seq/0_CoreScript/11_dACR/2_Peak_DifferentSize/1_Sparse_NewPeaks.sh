#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11.%j.err    # Standard error log
#SBATCH --array=0-12

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Combined_AnnV3_Onlychr.bed  \
  -b /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1/A619Bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_MergedPeak_Intergenic.bed  \
  -wa -wb \
  -sorted  | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl \
   - > /scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_PromoterIntergenic/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_PeaksCount_byA619Barcode.txt

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Combined_AnnV3_Onlychr.bed  \
 -b /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1/A619Bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_MergedPeak_Intergenic.bed  \
 -wa -wb \
 -sorted  | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl \
  - > /scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_PromoterIntergenic/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_PeaksCount_byBif3Barcode.txt
