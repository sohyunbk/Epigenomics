#!/bin/bash
#SBATCH --job-name=CPMNormalizationDAP        # Job name
#SBATCH --partition=schmitz_hm_p             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=04:50:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/CPMNormalizationDAP.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/CPMNormalizationDAP.%j.err    # Standard error log

ml Anaconda3/2023.09-0

cd /scratch/sb14489/7.DAPorChIP/DAPseq_WUS/WUS1_fastq_Mapped

conda activate Jbrowse

ml SAMtools/1.16.1-GCC-11.3.0
ml BEDTools/2.30.0-GCC-12.2.0

samtools sort -@24 -o HB67_WUS1_B73v5_Q30.sorted.bam HB67_WUS1_B73v5_Q30.bam
samtools index HB67_WUS1_B73v5_Q30.sorted.bam

##
#[bam_sort_core] merging from 0 files and 24 in-memory blocks...
#[main_samview] fail to read the header from "HB67_WUS1_B73v5_Q30.Sorted.coverage.bedGraph".
#(standard_in) 1: syntax error
#HB67_WUS1_B73v5_Q30.CPM.coverage.bedGraph is not case-sensitive sorted at line 22377802.  Ple
#ase use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again.

bedtools genomecov -ibam HB67_WUS1_B73v5_Q30.sorted.bam -bg > HB67_WUS1_B73v5_Q30.coverage.bedGraph
#sort -k1,1 -k2,2n HB67_WUS1_B73v5_Q30.coverage.bedGraph -o HB67_WUS1_B73v5_Q30.Sorted.coverage.bedGraph
TOTAL_READS=$(samtools view -F 260 -c HB67_WUS1_B73v5_Q30.sorted.bam)
#TOTAL_READS=$(samtools view -F 260 -c HB67_WUS1_B73v5_Q30.Sorted.coverage.bedGraph)
echo $TOTAL_READS
SCALE_FACTOR=$(echo "scale=6; 1/(${TOTAL_READS}/1000000)" | bc)
awk -v scale=${SCALE_FACTOR} 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*scale}' HB67_WUS1_B73v5_Q30.coverage.bedGraph > HB67_WUS1_B73v5_Q30.CPM.coverage.bedGraph
sort -k1,1 -k2,2n HB67_WUS1_B73v5_Q30.CPM.coverage.bedGraph -o HB67_WUS1_B73v5_Q30.CPM.coverage.Sorted.bedGraph
bedGraphToBigWig HB67_WUS1_B73v5_Q30.CPM.coverage.Sorted.bedGraph /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai HB67_WUS1_B73v5_Q30.CPM.bw
