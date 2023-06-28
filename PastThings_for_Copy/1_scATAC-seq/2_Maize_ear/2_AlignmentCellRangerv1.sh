#for i in 1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2
#for i in 2_rel2
for i in 2_rel2_2 3_bif3_2
do
echo $i
sbatch /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/2_Alignment/2-1_Alignment_CellRangerv1.sh /scratch/sb14489/3.scATAC/2.Maize_ear/ Maize_B73/Maize_B73_V5_withMtPt_scATACCellRangerv1  $i 2.Mapped_CellRangerv1
done
