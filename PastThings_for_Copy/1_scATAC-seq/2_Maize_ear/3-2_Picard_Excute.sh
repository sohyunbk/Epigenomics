for i in 1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2
#for i in 3_bif3 #Try new version
do
echo $i
sbatch /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/3_RefiningReads/3-2_PICARD_Cellrangerv1.sh /scratch/sb14489/3.scATAC/2.Maize_ear/  $i
done
