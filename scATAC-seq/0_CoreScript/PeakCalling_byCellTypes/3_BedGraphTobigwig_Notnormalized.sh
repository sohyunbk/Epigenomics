#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9-3_Bigwig.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9-3_Bigwig.%j.err    # Standard error log
#SBATCH --array=0-12

#f"bedGraphToBigWig {bed_graph_file_normalized} {reference_file} {bw_file_name}"
#bedGraphToBigWig Test_PeakCalling/Cluster2/Test_Cluster2.pool.macs_treat_pileup.normalized.bdg /scratch/sb14489/3.scATAC_flo/0.Reference/genome.fa.fai Test_PeakCalling/Cluster2/Test_Cluster2.normalized.bw

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc

#ClusterN=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
#ClusterN=(Bundle_Sheath Meta_proto_xylem Companion_cell Parenchyma Determinate_later_organs Phloem_sieve_element IM_SM_SPM Procambial_meristem L1_Layer Sclerenchyma Meta_proto_phloem_sieve_element SM_SPM_Base)
ClusterN=(Bundle_Sheath G2_M L1_Layer_of_Determinate_later_organs Parenchyma SM_SPM_Base Companion_cells_phloem_parenchyma_cells IM_SM_SPM Meta_proto_phloem_sieve_element Procambial_meristem Determinate_later_organs L1_Layer Meta_proto_xylem Sclerenchyma)

cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V2
mkdir BwFiles

# 1) A619
SampleName=A619
cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V2/"$SampleName"
bedSort  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.bdg  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup_Sorted.bdg
bedGraphToBigWig ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ../BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw

# 2) bif3
#SampleName=bif3
#cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V2/"$SampleName"
#bedSort  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.bdg  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup_Sorted.bdg
#bedGraphToBigWig ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ../BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw

# 3) bif3
#SampleName=rel2
#cd /scratch/sb14489/3.scATAC_flo/9_PeakCalling/"$SampleName"_PeakCalling
#bedSort ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized.bdg  ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized_2Sorted.bdg
#bedGraphToBigWig ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized_2Sorted.bdg /scratch/sb14489/3.scATAC_flo/0.Reference/genome.fa.fai ../BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".normalized.bw
#bedGraphToBigWig Test.bdg  /scratch/sb14489/3.scATAC_flo/0.Reference/genome.fa.fai Test.bw

# 4) bif3
#SampleName=relk1
#cd /scratch/sb14489/3.scATAC_flo/9_PeakCalling/"$SampleName"_PeakCalling
#bedSort ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized.bdg ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized_2Sorted.bdg
#bedGraphToBigWig ./"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".pool.macs_treat_pileup.normalized_2Sorted.bdg /scratch/sb14489/3.scATAC_flo/0.Reference/genome.fa.fai ../BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".normalized.bw
