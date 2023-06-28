args = get_parser().parse_args()
GenomeSize =  "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"
GenomeSize_Dic = Make_Fai_Dic(GenomeSize)

Dir = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619" ## Inputfiles dir

SampleName = "A619"
PeakLength = int(200)

OutputDir = "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1"

python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8-2_Peaks_Cut/MergeCutPeaks_toFixHumpedPeaks.py -method Method1 \
    -inputpath1 /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619 \
     -inputpath2 /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3 \
      -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
      -OutputFileName A619Bif3 -OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1 \
      -Ann /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed \
      -Sparse1 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Combined.bed \
      -Sparse2 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Combined.bed


python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8-2_Peaks_Cut/MergeCutPeaks_toFixHumpedPeaks.py -method Method1 \
 -inputpath1 /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619 \
 -inputpath2 /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3 \
 -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
  -OutputFileName A619Bif3 -OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1 \
  -Ann /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed \
  -Sparse1 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/A619_Ex.sparse \
  -Sparse2 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/Bif3_Ex.sparse



python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8-2_Peaks_Cut/MergeCutPeaks_toFixHumpedPeaks.py -method Method2 \
    -inputpath /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619 \
    -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
    -OutputFileName A619 -OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Method2 \
    -size 200
