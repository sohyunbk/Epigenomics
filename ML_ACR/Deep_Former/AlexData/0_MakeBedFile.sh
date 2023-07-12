cd /scratch/sb14489/8.ML_ACR/1.InputBed

sort -k1V -k2n -k3n Seedling_Peaks.bed > sorted_Seedling_Peaks.bed

module load tabix/0.2.6-GCCcore-8.3.0

bgzip -c sorted_Seedling_Peaks.bed > sorted_Seedling_Peaks.bed.gz
tabix -p bed sorted_Seedling_Peaks.bed.gz
cut -f 4 sorted_Seedling_Peaks.bed | sort -u > Seedling_distinct_features.txt

bgzip -c sorted_Seedling_Peaks_200bp.bed > sorted_Seedling_Peaks_200bp.bed.gz
tabix -p bed sorted_Seedling_Peaks_200bp.bed.gz

bgzip -c sorted_Seedling_Peaks_204bp.bed > sorted_Seedling_Peaks_204bp.bed.gz
tabix -p bed sorted_Seedling_Peaks_204bp.bed.gz

### New bed files
sort -k1V -k2n -k3n  /scratch/sb14489/8.ML_ACR/1.InputBed/AllACRs_130539.500bp.bed > /scratch/sb14489/8.ML_ACR/1.InputBed/AllACRs_130539.500bp.Sorted.bed
sort -k1V -k2n -k3n  /scratch/sb14489/8.ML_ACR/1.InputBed/celltype_RestrictACRs_118160.500bp.bed > /scratch/sb14489/8.ML_ACR/1.InputBed/celltype_RestrictACRs_118160.500bp.Sorted.bed

bgzip -c /scratch/sb14489/8.ML_ACR/1.InputBed/AllACRs_130539.500bp.Sorted.bed > /scratch/sb14489/8.ML_ACR/1.InputBed/AllACRs_130539.500bp.Sorted.bed.gz
tabix -p bed /scratch/sb14489/8.ML_ACR/1.InputBed/AllACRs_130539.500bp.Sorted.bed.gz

bgzip -c /scratch/sb14489/8.ML_ACR/1.InputBed/celltype_RestrictACRs_118160.500bp.Sorted.bed >  /scratch/sb14489/8.ML_ACR/1.InputBed/celltype_RestrictACRs_118160.500bp.Sorted.bed.gz
tabix -p bed /scratch/sb14489/8.ML_ACR/1.InputBed/celltype_RestrictACRs_118160.500bp.Sorted.bed.gz
