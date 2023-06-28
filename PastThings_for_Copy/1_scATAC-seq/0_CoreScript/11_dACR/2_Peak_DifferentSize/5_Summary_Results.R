## Summarize the results and make bed file for 

#### Files Changed ###
ResultFileDir <- "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/"
NameRules <- ".EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
#################
CellType <- c(
"BundleSheath_VascularSchrenchyma",
"CalloseRelated",
"FloralMeristem_SuppressedBract",
"G2_M", 
"IM-OC",
"IM_SPM_SM",
"L1",
"L1atFloralMeristem",
"PhloemPrecursor",
"ProcambialMeristem_ProtoXylem_MetaXylem",
"ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma",
"SPM-base_SM-base",
"XylemParenchyma_PithParenchyma")

CT <- "IM-OC"

Summary <- data.frame()

for (CT in CellType){
  #EdgeRResult <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/",
  #                                       CT,".EdgeRResult_PseudoReplicate_withPromoterRegion.txt"),header=TRUE)
  EdgeRResult <- read.table(paste0(ResultFileDir,
                                   CT,NameRules),header=TRUE)
  Total <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3_MergePeakbyCelltypes_Method1/A619Bif3_",
                             CT ,"_MergedPeak_Intergenic.bed"))
  #FisherResult <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks/",CT,".FisherExactTest.txt"),header=TRUE)
  
  head(EdgeRResult)
  Sig <- EdgeRResult[EdgeRResult$FDR < 0.05,]
  #Sig <- EdgeRResult[EdgeRResult_Mixed$FDR < 0.05,]
  
  nrow(Sig)
  tail(Sig)
  LogFC_positive <- nrow(Sig[Sig$logFC > 0,])
  LogFC_negative <- nrow(Sig[Sig$logFC < 0,])
  #length(LogFC[LogFC== TRUE])
  #length(LogFC[LogFC== FALSE])
  Sig
  #head(FisherResult)
  
  #Sig2 <- FisherResult[FisherResult$FDR < 0.05,]
  #nrow(Sig2)
  #nrow(Total)
  #temp <- data.frame(nrow(Total),nrow(Sig), nrow(Sig2))
  temp <- data.frame(CT,nrow(Total),nrow(Sig),LogFC_positive,LogFC_negative)
  Summary <- rbind(Summary,temp)
  ## positive logFC means Bif3 is higher
  ## negative logFC means A619 is higher
}
write.table(Summary,paste0(ResultFileDir,"/Summary.txt"),sep="\t", quote=F, row.names=F)
