module load BEDTools/2.30.0-GCC-10.2.0
bedtools getfasta -fi /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa -bed IM-OC.A619Higher.Bed > IM-OC.A619Higher.fa
bedtools getfasta -fi /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa -bed IM-OC.Bif3Higher.Bed > IM-OC.Bif3Higher.fa


module load MEME/5.4.1-foss-2019b-Python-3.7.4
meme-chip IM-OC.Bif3Higher.fa -minw 4 -maxw 15 -o IM-OC.Bif3Higher
