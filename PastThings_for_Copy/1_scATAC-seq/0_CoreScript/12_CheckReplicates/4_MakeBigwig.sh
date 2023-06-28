#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Bigwig.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Bigwig.%j.err    # Standard error log
#SBATCH --array=0-12

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

cd /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates
mkdir BwFiles

SampleName=3_bif3_Re1
bedSort  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM.bdg  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg
bedGraphToBigWig ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ./BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw

SampleName=3_bif3_Re2
bedSort  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM.bdg  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg
bedGraphToBigWig ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ./BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw

SampleName=1_A619_Re1
bedSort  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM.bdg  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg
bedGraphToBigWig ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ./BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw

SampleName=1_A619_Re2
bedSort  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM.bdg  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg
bedGraphToBigWig ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ./BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw
