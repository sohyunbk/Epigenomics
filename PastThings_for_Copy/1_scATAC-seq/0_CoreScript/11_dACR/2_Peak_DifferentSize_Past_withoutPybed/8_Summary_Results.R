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
  EdgeRResult <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks/",CT,".EdgeRResult_PseudoReplicate.txt"),header=TRUE)
  Total <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_For_dACR/",CT,"_ComA619Bif3_BLRemove_Intergenic.bed"))
  FisherResult <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks/",CT,".FisherExactTest.txt"),header=TRUE)
  
  head(EdgeRResult)
  Sig <- EdgeRResult[EdgeRResult$FDR < 0.05,]
  nrow(Sig)
  tail(Sig)
  LogFC <- (Sig$logFC > 0)
  length(LogFC[LogFC== TRUE])
  length(LogFC[LogFC== FALSE])
  
  head(FisherResult)
  Sig2 <- FisherResult[FisherResult$FDR < 0.05,]
  nrow(Sig2)
  nrow(Total)
  temp <- data.frame(nrow(Total),nrow(Sig), nrow(Sig2))
  Summary <- rbind(Summary,temp)
  ## positive logFC means Bif3 is higher
  ## negative logFC means A619 is higher
  }
