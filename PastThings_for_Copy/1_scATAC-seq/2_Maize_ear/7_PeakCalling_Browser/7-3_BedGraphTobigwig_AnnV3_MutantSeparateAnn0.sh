#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9-3_Bigwig.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9-3_Bigwig.%j.err    # Standard error log
#SBATCH --array=0-12

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc

ClusterN=(BundleSheath_VascularSchrenchyma L1atFloralMeristem CalloseRelated PhloemPrecursor FloralMeristem_SuppressedBract ProcambialMeristem_ProtoXylem_MetaXylem G2_M ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma IM-OC SPM-base_SM-base IM_SPM_SM XylemParenchyma_PithParenchyma L1)

AnnDir=Ann_V3_SeparateAnn_Mutant
cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/"$AnnDir"/
mkdir BwFiles

# bif3
SampleName=bif3
cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/"$AnnDir"/"$SampleName"
bedSort  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized.bdg  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized_2Sorted.bdg
bedGraphToBigWig ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized_2Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ../BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".normalized.bw
