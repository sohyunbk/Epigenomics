#!/bin/bash
#SBATCH --job-name=11_ChangeFileForJbrowse        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=30:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/12_Macs2_replicates.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/12_Macs2_replicates.%j.err    # Standard error log
#SBATCH --array=0-12                  # Array range

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

macs2 callpeak -t /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/1_A619_Re1_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bed -f BED --nomodel \
      --keep-dup all --extsize 150 --shift -50 --qvalue .05 \
      --outdir /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/1_A619_Re1_"${ClusterN[SLURM_ARRAY_TASK_ID]}"  \
      --bdg -n 1_A619_Re1_"${ClusterN[SLURM_ARRAY_TASK_ID]}"


macs2 callpeak -t /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/1_A619_Re2_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bed -f BED --nomodel \
      --keep-dup all --extsize 150 --shift -50 --qvalue .05 \
      --outdir /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/1_A619_Re2_"${ClusterN[SLURM_ARRAY_TASK_ID]}"  \
      --bdg -n 1_A619_Re2_"${ClusterN[SLURM_ARRAY_TASK_ID]}"

macs2 callpeak -t /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/3_bif3_Re1_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bed -f BED --nomodel \
      --keep-dup all --extsize 150 --shift -50 --qvalue .05 \
      --outdir /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/3_bif3_Re1_"${ClusterN[SLURM_ARRAY_TASK_ID]}"  \
      --bdg -n 3_bif3_Re1_"${ClusterN[SLURM_ARRAY_TASK_ID]}"


macs2 callpeak -t /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/3_bif3_Re2_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bed -f BED --nomodel \
      --keep-dup all --extsize 150 --shift -50 --qvalue .05 \
      --outdir /scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/3_bif3_Re2_"${ClusterN[SLURM_ARRAY_TASK_ID]}"  \
      --bdg -n 3_bif3_Re2_"${ClusterN[SLURM_ARRAY_TASK_ID]}"
