#!/bin/bash
#SBATCH --job-name=MetaPlot        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=01:50:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/MetaPlot.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MetaPlot.%j.err    # Standard error log


cd /scratch/sb14489/7.DAPorChIP/DAPseq_WUS/WUS1_fastq_Mapped

module load deepTools/3.5.2-foss-2022a

computeMatrix reference-point -S ./HB67_WUS1_B73v5_Q30_bl.bigwig \
                            -R /scratch/sb14489/7.DAPorChIP/DAPseq_WUS/HB67_WUS1_B73v5_Q30_qval5_finalBl/HB67_WUS1_B73v5_Q30_qval5_finalBl.GEM_events.narrowPeak \
                              --beforeRegionStartLength 3000 \
                              --afterRegionStartLength 3000 \
                              --referencePoint center \
                               -o ./ZmWUS1.mat.gz \
				                       --binSize 10 \
                               -p 20


#plotProfile -m ./7.MetaPlot/SoybeanARF"${List[SLURM_ARRAY_TASK_ID]}".mat.gz \
#              -out ./7.MetaPlot/SoybeanARFs"${List[SLURM_ARRAY_TASK_ID]}".png \
#             --colors "#cf7c4c" \
#              --refPointLabel "Peak Summit" \
#               --plotTitle "" \
#               --samplesLabel "ARF${List[SLURM_ARRAY_TASK_ID]}" \
#               --regionsLabel "Read"
