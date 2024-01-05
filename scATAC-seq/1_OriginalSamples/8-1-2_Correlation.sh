#!/bin/bash
#SBATCH --job-name=Correlation        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Correlation.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Correlation.%j.err    # Standard error log
#SBATCH --array=0-1


source activate r_env
S1SampleSparse=(
/A619_Bif3_500bpCommonPeak/ComA619Bif3_A619Barcode_Tn5Count.sparse
/A619_rel2_500bpCommonPeak/ComPeakA619rel2_A619Barcode_Tn5Count.sparse)
S2SampleSparse=(
/A619_Bif3_500bpCommonPeak/ComA619Bif3_Bif3Barcode_Tn5Count.sparse
/A619_rel2_500bpCommonPeak/ComPeakA619rel2_Rel2Barcode_Tn5Count.sparse)
S2SampleName=(
Bif3
rel2
)
S2Metas=(
/Bif3/Bif3_AnnV4_metadata.txt
/rel2/rel2_AnnV4.txt
)
CommonPeak=(
/A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed
/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed
)
CommonPeakIntergenic=(
/A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks_Intergenic.bed
A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks_Intergenic.bed
)
OutfileName=(
A619andBif3
A619andrel2
)
CTNameOrder=(
/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt
/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forrel2.txt
)
Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Correlation/Correlation_Intergenic2000MostVariationACR_500pbCommonACR.R \
  --S1_Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${S1SampleSparse[SLURM_ARRAY_TASK_ID]}" \
  --S2_Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${S2SampleSparse[SLURM_ARRAY_TASK_ID]}" \
  --S1Name A619 \
  --S2Name "${S2SampleName[SLURM_ARRAY_TASK_ID]}" \
  --S1_Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt \
  --S2_Meta //scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${S2Metas[SLURM_ARRAY_TASK_ID]}" \
  --ClusterColumnName Ann_v4 \
  --S1and2_500bpPeak //scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${CommonPeak[SLURM_ARRAY_TASK_ID]}" \
  --S1and2_500bpInterPeak //scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${CommonPeakIntergenic[SLURM_ARRAY_TASK_ID]}" \
  --OutPath /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/1.Correlation \
  --OutFileName "${OutfileName[SLURM_ARRAY_TASK_ID]}"
