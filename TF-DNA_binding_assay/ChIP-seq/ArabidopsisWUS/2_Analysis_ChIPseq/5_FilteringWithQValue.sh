# 8th Column: pValue -log10
# 9th Column: QValue -log10
#!/bin/bash

cd /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling_HighCutoff
input_file="WUS_GS_Q0.001_peaks.narrowPeak"

for threshold in 10 20 30 40 50
do
    output_file_Chr="WUS_GS_Q${threshold}.narrowPeak_Chr"
    output_file="WUS_GS_Q${threshold}.narrowPeak"
    # Use awk to filter lines where the 9th column is greater than the current threshold
    awk -F'\t' -v threshold="$threshold" '$9 > threshold {print "Chr"$0}' "$input_file" > "$output_file_Chr"
    awk -F'\t' -v threshold="$threshold" '$9 > threshold {print $0}' "$input_file" > "$output_file"
done
