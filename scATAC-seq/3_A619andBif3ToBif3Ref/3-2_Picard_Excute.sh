for i in 1_A619 3_bif3 1_A619_2 3_bif3_2
#for i in 3_bif3 #Try new version
do
echo $i
sbatch /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/3_RefiningReads/3-2_PICARD.sh /scratch/sb14489/3.scATAC/4.Bif3Ref/  $i
done
