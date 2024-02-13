#!/bin/bash
#SBATCH --job-name=ChromVar        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ChromVar.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ChromVar.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-3

ml Anaconda3/2020.02
source activate JASPAR_act

WDs=(
AnnV4_CommonA619Bif3Peak
AnnV4_CommonA619Bif3Peak
AnnV4
AnnV4
)

Sparses=(
/A619_Bif3_500bpCommonPeak/ComA619Bif3_A619Barcode_Tn5Count.sparse
/A619_Bif3_500bpCommonPeak/ComA619Bif3_A619Barcode_Tn5Count.sparse
/A619/A619Intergenic500bp_A619Barcode_Tn5Count.sparse
Bif3/Bif3Intergenic500bp_Bif3Barcode_Tn5Count.sparse
)

Metas=(
A619/Ref_AnnV4_metadata.txt
Bif3/Bif3_AnnV4_metadata.txt
A619/Ref_AnnV4_metadata.txt
Bif3/Bif3_AnnV4_metadata.txt
)

PCSs=(
Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt
Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt
Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt
)
Markovs=(
/A619_IncludingZmCLE7/A619_IncludingZmCLE7.MarkovMatrix.rds
/Bif3_IncludingZmCLE7/Bif3_IncludingZmCLE7.MarkovMatrix.rds
/A619_IncludingZmCLE7/A619_IncludingZmCLE7.MarkovMatrix.rds
/Bif3_IncludingZmCLE7/Bif3_IncludingZmCLE7.MarkovMatrix.rds
)
Samples=(
A619
Bif3
A619
Bif3
)
IntergenicPeaks=(
A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks_Intergenic.bed
A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks_Intergenic.bed
A619/A619.500bp_peaks_Intergenic.bed
Bif3/Bif3.500bp_peaks_Intergenic.bed
)
Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MotifDeviation/MotifAnalysis_ChromVar_AllMotifs_SmoothData.R  \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/"${WDs[SLURM_ARRAY_TASK_ID]}"  \
 --Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${Sparses[SLURM_ARRAY_TASK_ID]}" \
 --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${Metas[SLURM_ARRAY_TASK_ID]}" \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/"${PCSs[SLURM_ARRAY_TASK_ID]}" \
 --Markov /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/"${Markovs[SLURM_ARRAY_TASK_ID]}" \
 --SampleName A619 --IGPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${IntergenicPeaks[SLURM_ARRAY_TASK_ID]}"
