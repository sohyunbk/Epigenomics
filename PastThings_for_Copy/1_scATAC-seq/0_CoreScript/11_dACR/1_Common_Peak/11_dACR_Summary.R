setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_OnlyIntergenic")

Files <- list.files(pattern=".EdgeRResult.txt")
for (sFile in Files){
  print(sFile)
  Table <- read.table(sFile, header=TRUE)
  print(sum(Table$FDR < 0.05))
}
Table <- read.table("IM-OC.EdgeRResult.txt", header=TRUE)
head(Table)
Sig <- Table[Table$FDR < 0.05,]
dim(Sig)
sum(Sig$logFC > 0)
sum(Sig$logFC < 0)
head(Sig)
Sig[Sig$logFC > 0,]

head(Sig)
library(tidyr)
Bedfile_Under0.5FDR <- separate(Sig,col=Peak,into=c("Chr","Start","End"),sep="_")
write.table(Bedfile_Under0.5FDR, file="IM-OC.FDR0.5.bed", quote=F, row.names=F, col.names=F, sep="\t")

#bedtools intersect -wo -a IM-OC.FDR0.5.bed -b ../../../../7.DAPorChIP/DAPseq_WUS/HB67_WUS1_B73v5_Q30_qval5_finalBl/HB67_WUS1_B73v5_Q30_qval5_finalBl.GEM_events.narrowPeak > Overlap_IM-OC.FDR0.5_WUS1Q5.txt