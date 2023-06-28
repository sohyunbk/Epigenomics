#!/bin/bash
#SBATCH --job-name=11_ChangeFileForJbrowse        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11_ChangeFileForJbrowse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11_ChangeFileForJbrowse.%j.err    # Standard error log
#SBATCH --array=0-12


cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/Jbrowse

SampleName="A619"
/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed ./"$SampleName"/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks --trackLabel "$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks --out ./
/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed ./"$SampleName"/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_summits --trackLabel "$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_summits --out ./


SampleName="bif3"
/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed ./"$SampleName"/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks --trackLabel "$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks --out ./
/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed ./"$SampleName"/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_summits --trackLabel "$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_summits --out ./

#scp -r bif3_IM-OC.reproducible_narrow_peaks schmitzlab1@heredity.genetics.uga.edu:/data01/epigenome/JBrowse/maize_v5_UpdatedAnn/tracks
#Schmacct5$
