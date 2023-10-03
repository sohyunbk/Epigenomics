cd /scratch/sb14489/8.ML_ACR/2.MaizeEar/1.InputBed

sort -k1V -k2n -k3n MaizeEar_500bpPeak_A619_byCT.bed > sorted_MaizeEar_500bpPeak_A619_byCT.bed

module load tabix/0.2.6-GCCcore-11.3.0

bgzip -c sorted_MaizeEar_500bpPeak_A619_byCT.bed > sorted_MaizeEar_500bpPeak_A619_byCT.bed.gz
tabix -p bed sorted_MaizeEar_500bpPeak_A619_byCT.bed.gz
cut -f 4 sorted_MaizeEar_500bpPeak_A619_byCT.bed | sort -u > Ear_A619_features.txt
