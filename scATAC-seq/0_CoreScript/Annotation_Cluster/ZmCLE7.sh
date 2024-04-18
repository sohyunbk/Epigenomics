cd /scratch/sb14489/0.Reference/Maize_B73/v3
module load BEDTools/2.30.0-GCC-12.2.0
module load Cufflinks/20190706-GCC-11.2.0
gffread -w Zea_mays.AGPv3.19_onlyZmCLE7_Promoter.fa \
-g B73_RefGen_v3.fa Zea_mays.AGPv3.19_onlyZmCLE7.gtf

bedtools getfasta -fi B73_RefGen_v3.fa -bed Zea_mays.AGPv3.19_onlyZmCLE7.gtf \
 -fo Zea_mays.AGPv3.19_onlyZmCLE7_Promoter.fa

 bedtools getfasta -fi B73_RefGen_v3.fa -bed Zea_mays.AGPv3.19_onlyZmCLE7.bed \
  -fo Zea_mays.AGPv3.19_onlyZmCLE7_Promoter_Bed.fa


 bedtools getfasta -fi ../Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa -bed Test_ZmCLE7.bed    \
  -fo Test_ZmCLE7.fa
