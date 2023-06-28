args <- commandArgs(T)
CT <- as.character(args[1])

Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks/"
CountTable <- readRDS(paste0(Path,CT,"_CountMatrix_PerCelltype.rds"))
CellNames <- colnames(CountTable)

Result <- data.frame()
for (nNumb in c(1:length(CountTable$PeakLocus))){
  Peak <- CountTable$PeakLocus[nNumb]
  CountPart <- CountTable[nNumb,]
  #CountPart["CB:Z:AATGTCGGTCTCTGGG-3_bif3_2"]
  A619_CountPart <- CountPart[,endsWith(CellNames, 'A619')]
  A619_2_CountPart <- CountPart[,endsWith(CellNames, 'A619_2')]
  A619_Combined <- cbind(A619_CountPart,A619_2_CountPart)
  A619_table <- table(as.numeric(A619_Combined[1,]))
  NonTn5_A619 <- as.numeric(A619_table[1])
  Tn5_A619 <- length(as.numeric(A619_Combined[1,])) - NonTn5_A619
  
  bif3_CountPart <- CountPart[,endsWith(CellNames, 'bif3')]
  bif3_2_CountPart <- CountPart[,endsWith(CellNames, 'bif3_2')]
  bif3_Combined <- cbind(bif3_CountPart,bif3_2_CountPart)
  
  bif3_table <- table(as.numeric(bif3_Combined[1,]))
  NonTn5_bif3<- as.numeric(bif3_table[1])
  Tn5_bif3<- length(as.numeric(bif3_Combined[1,])) - NonTn5_bif3
  
  dat <- data.frame(
    "A619" = c(NonTn5_A619, Tn5_A619),
    "Bif3" = c(NonTn5_bif3, Tn5_bif3),
    row.names = c("NoneTn5", "Tn5"),
    stringsAsFactors = FALSE
  )
  colnames(dat) <- c("A619", "Bif3")
  test <- fisher.test(dat)
  TempTable <- data.frame(Peak,as.numeric(test$p.value),as.numeric(test$estimate),NonTn5_A619,Tn5_A619,NonTn5_bif3,Tn5_bif3)
  Result <- rbind(Result,TempTable)

}

colnames(Result) <- c("Peak","P.Value","OddRatio","NonTn5_A619","Tn5_A619","NonTn5_bif3","Tn5_bif3")
Result$FDR <- p.adjust(Result$P.Value,method="hochberg")

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks")
write.table(Result, file=paste0(CT,".FisherExactTest.txt"), quote=F, row.names=F, col.names=T, sep="\t")
