#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11_intersect.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11_intersect.%j.err    # Standard error log
#SBATCH --array=0-12

ClusterN=(BundleSheath_VascularSchrenchyma CalloseRelated FloralMeristem_SuppressedBract G2_M IM-OC IM_SPM_SM L1 L1atFloralMeristem PhloemPrecursor ProcambialMeristem_ProtoXylem_MetaXylem ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma SPM-base_SM-base XylemParenchyma_PithParenchyma)

cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak

module load BEDTools/2.30.0-GCC-10.2.0

## 1) Grep only Chr : as bedtools makes error if chr - Mt not matched
grep  '^chr' ./A619/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/A619_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks >  ./A619/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/A619_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks_Onlychr
grep  '^chr' ./bif3/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks >  ./bif3/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks_Onlychr

## 2) intersect
bedtools intersect -wo -a  ./A619/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/A619_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks_Onlychr \
  -b ./bif3/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/bif3_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks_Onlychr  > \
 ./A619_bif3_For_dACR/"${ClusterN[SLURM_ARRAY_TASK_ID]}"_A619_bif3_Intersect.txt
