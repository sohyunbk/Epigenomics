#!/bin/bash
#SBATCH --job-name=ChangeFileForJbrowse        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ChangeFileForJbrowse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ChangeFileForJbrowse.%j.err    # Standard error log
#SBATCH --array=0-12

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/Jbrowse

### All ACR A619_list
cd /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1
/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed A619Bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_MergedPeak.bed \
 --trackLabel "${ClusterN[SLURM_ARRAY_TASK_ID]}".A619Bif3_MergedPeakByCT --out ./

### dACR list
cd /scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/
cat "${ClusterN[SLURM_ARRAY_TASK_ID]}".A619Higher.Bed "${ClusterN[SLURM_ARRAY_TASK_ID]}".Bif3Higher.Bed > "${ClusterN[SLURM_ARRAY_TASK_ID]}".FDR0.05.Bed
/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed "${ClusterN[SLURM_ARRAY_TASK_ID]}".FDR0.05.Bed \
 --trackLabel "${ClusterN[SLURM_ARRAY_TASK_ID]}".FDR0.05 --out ./
