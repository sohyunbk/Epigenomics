#!/bin/bash
#SBATCH --job-name=MetaPlot        # Job name
#SBATCH --partition=schmitz_hm_p             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=01:50:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/MetaPlot.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MetaPlot.%j.err    # Standard error log

ml SAMtools/1.16.1-GCC-11.3.0
ml BEDTools/2.30.0-GCC-12.2.0

cd /scratch/sb14489/7.DAPorChIP/DAPseq_WUS/WUS1_fastq_Mapped

source  acitivate Jbrowse

samtools sort -@24 -o HB67_WUS1_B73v5_Q30.sorted.bam HB67_WUS1_B73v5_Q30.bam
samtools index HB67_WUS1_B73v5_Q30.sorted.bam

##
bedtools genomecov -ibam HB67_WUS1_B73v5_Q30.sorted.bam -bg > HB67_WUS1_B73v5_Q30.coverage.bedGraph


TOTAL_READS=$(samtools view -F 260 -c HB67_WUS1_B73v5_Q30.sorted.bam)
echo $TOTAL_READS
SCALE_FACTOR=$(echo "scale=6; 1/(${TOTAL_READS}/1000000)" | bc)
awk -v scale=${SCALE_FACTOR} 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*scale}' HB67_WUS1_B73v5_Q30.coverage.bedGraph > HB67_WUS1_B73v5_Q30.CPM.coverage.bedGraph
bedGraphToBigWig HB67_WUS1_B73v5_Q30.CPM.coverage.bedGraph chrom.sizes HB67_WUS1_B73v5_Q30.CPM.bw
