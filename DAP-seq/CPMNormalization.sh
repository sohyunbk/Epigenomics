ml SAMtools/1.16.1-GCC-11.3.0
samtools sort -@4 -o HB67_WUS1_B73v5_Q30.sorted.bam HB67_WUS1_B73v5_Q30.bam
samtools index HB67_WUS1_B73v5_Q30.sorted.bam
genomeCoverageBed -ibam HB67_WUS1_B73v5_Q30.sorted.bam -bg -split > coverage.bedGraph
SCALE_FACTOR=$(echo "scale=6; 1/(${TOTAL_MAPPED_READS}/1000000)" | bc)
awk -v scale=${SCALE_FACTOR} 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*scale}' coverage.bedGraph > scaled_coverage.bedGraph
bedGraphToBigWig scaled_coverage.bedGraph chrom.sizes output.bw
