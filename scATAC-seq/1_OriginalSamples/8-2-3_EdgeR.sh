#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11.%j.err    # Standard error log
#SBATCH --array=0-13

## 14 number of cell types
ClusterN=(L1atFloralMeristem  Unknown1 FloralMeristem_SuppressedBract  PhloemPrecursor  Unknown2
G2_M ProcambialMeristem_ProtoXylem_MetaXylem  Unknown_Sclerenchyma IM-OC ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
Unknown_lowFRiP L1 SPM-base_SM-base  XylemParenchyma_PithParenchyma)

bedFiles=(1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr
1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr
1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr
1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr
1_A619_Combined_Sorted_k12.bed_OnlyChr 1_A619_Combined_Sorted_k12.bed_OnlyChr
3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr
3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr
 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr
 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr
 3_bif3_Combined_Sorted_k12.bed_OnlyChr 3_bif3_Combined_Sorted_k12.bed_OnlyChr
)
OutFileName=(byA619Barcode byA619Barcode byA619Barcode byA619Barcode byA619Barcode byA619Barcode byA619Barcode byA619Barcode
byA619Barcode byA619Barcode byA619Barcode byA619Barcode byA619Barcode byA619Barcode
byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode
byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode byBif3Barcode)

ml Anaconda3/2023.09-0
source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Correlation/Correlation_Intergenic2000MostVariationACR_500pbCommonACR.R \
