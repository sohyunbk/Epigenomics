##### Back to EdgeR!! 
#conda activate r4-base
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)

## 1. Load files and get filtered Sparse file for save space
## --> should combine A619+Bif3 in the begining to keep all the peaks as features

#CT <- "XylemParenchyma_PithParenchyma"
args <- commandArgs(T)
CT <- as.character(args[1])
CT <-"IM-OC"

Sparsefile_A619 <- paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_For_dACR/",CT ,"_A619_toPeaks.sparse")
Sparsefile_Bif3  <- paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_For_dACR/",CT ,"_Bif3_toPeaks.sparse")
MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
InterGenic_peak <- read.table(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_For_dACR/",CT ,"_ComA619Bif3_BLRemove_TSSEnd.bed"))
#chr10   149311748       149312391
dim(InterGenic_peak)
InterGenic_peak_Pos <- paste(InterGenic_peak$V1,InterGenic_peak$V2,InterGenic_peak$V3,sep="_")
head(InterGenic_peak_Pos)

GetFilteredSparseData <- function(Sparsefile,MetaFile,SelectedPeaksPos,CellType){
Sparse<- read.table(Sparsefile,header=F)
meta_data <- read.table(MetaFile, header=TRUE)
#head(meta_data)
meta_data <- meta_data[meta_data$Ann_v3 == CellType,]
#dim(meta_data)
colnames(Sparse) <- c("PeakLocus","cellID","accessability")
Sparse_Intergenic <- Sparse[Sparse$PeakLocus%in%SelectedPeaksPos,]
Sparse_Intergenic_SelectedCells <-Sparse_Intergenic[Sparse_Intergenic$cellID %in% meta_data$cellID,]
dim(Sparse_Intergenic_SelectedCells)
return(Sparse_Intergenic_SelectedCells)
}

Sparse_A619 <- GetFilteredSparseData(Sparsefile_A619,MetaFileA619,InterGenic_peak_Pos,CT)
Sparse_Bif3 <- GetFilteredSparseData(Sparsefile_Bif3,MetaFileBif3,InterGenic_peak_Pos,CT)

Sparse_Combined <- rbind(Sparse_A619,Sparse_Bif3)
Peak_Cell_Count <-  spread(Sparse_Combined,key = cellID,value =accessability)
dim(Peak_Cell_Count)
Peak_Cell_Count[is.na(Peak_Cell_Count)] <- 0

head(Peak_Cell_Count)[,c(1:10)]
#Peak_Cell_Count <- readRDS("Peak_perCell_Counts.rds")
Pos <- Peak_Cell_Count$PeakLocus

meta_data_A619 <- read.table(MetaFileA619, header=TRUE)
meta_data_Bif3 <- read.table(MetaFileBif3, header=TRUE)
CombinedMeta <- rbind(meta_data_A619,meta_data_Bif3)
head(CombinedMeta)

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks")

saveRDS(Peak_Cell_Count, file=paste0(CT,"_CountMatrix_PerCelltype.rds"))
Peak_Cell_Count <- readRDS(paste0(CT,"_CountMatrix_PerCelltype.rds"))

### 2. EdgeR

#### Make pseudobulk
DevideSampleIntoTwo <- function(StringVector){
  nSample <- length(StringVector)
  RandomSampleNumber <- sample(1:nSample, round(nSample/2), replace=F)
  Re1 <- StringVector[RandomSampleNumber]
  Re2 <- StringVector[!StringVector%in%Re1]
  return(list(Re1,Re2))
}
  
# Start with first cell type
Count_Celltype <- Peak_Cell_Count
CellNames <- colnames(Count_Celltype)

A619 <- CellNames[endsWith(CellNames, 'A619')]
A619_2 <- CellNames[endsWith(CellNames, 'A619_2')]
Bif3 <- CellNames[endsWith(CellNames, 'bif3')]
Bif3_2 <- CellNames[endsWith(CellNames, 'bif3_2')]

A619_list <- DevideSampleIntoTwo(A619)
A619_1_Re1 <- A619_list[[1]]
A619_1_Re2 <- A619_list[[2]]
A619_2_list <- DevideSampleIntoTwo(A619_2)
A619_2_Re1 <- A619_2_list[[1]]
A619_2_Re2 <- A619_2_list[[2]]

Bif3_list <- DevideSampleIntoTwo(Bif3)
Bif3_1_Re1 <- Bif3_list[[1]]
Bif3_1_Re2 <- Bif3_list[[2]]
Bif3_2_list <- DevideSampleIntoTwo(Bif3_2)
Bif3_2_Re1 <- Bif3_2_list[[1]]
Bif3_2_Re2 <- Bif3_2_list[[2]]

ReSampleList <- list()
ReSampleList[["A619_1_Re1"]] <-A619_1_Re1
ReSampleList[["A619_1_Re2"]] <-A619_1_Re2
ReSampleList[["A619_2_Re1"]] <-A619_2_Re1
ReSampleList[["A619_2_Re2"]] <-A619_2_Re2  

ReSampleList[["Bif3_1_Re1"]] <-Bif3_1_Re1
ReSampleList[["Bif3_1_Re2"]] <-Bif3_1_Re2
ReSampleList[["Bif3_2_Re1"]] <-Bif3_2_Re1
ReSampleList[["Bif3_2_Re2"]] <-Bif3_2_Re2  

Count_sPseudoReplicates <- data.frame(Count_Celltype[,1])
for (sPseudo in names(ReSampleList)){
  Count <- Count_Celltype[,ReSampleList[[sPseudo]]]
  Count_sPseudoReplicates <- cbind(Count_sPseudoReplicates,rowSums(Count))
}
Count_sPseudoReplicates <- Count_sPseudoReplicates[,-1]
colnames(Count_sPseudoReplicates) <-names(ReSampleList)
head(Count_sPseudoReplicates)

y <- DGEList(counts=Count_sPseudoReplicates, gene=Pos)
y <- calcNormFactors(y)

SampleName <- c("A619","A619","A619","A619","Bif3","Bif3","Bif3","Bif3")

design <- model.matrix(~SampleName)
y <- estimateDisp(y,design) ## It takes long time... 
fit <- glmFit(y, design) ## It takes long time...too.... # likelihood ratio test

TMM <- cpm(y, normalized.lib.sizes=TRUE,log=T)
TMMcolName <- colnames(TMM)
TMMcolName <- paste("TMM_",TMMcolName,sep="")

#str(TMMcolName)
colnames(TMM) <- as.factor(TMMcolName)

lrt <- glmLRT(fit,coef=2)
Result <- topTags(lrt, n=dim(Count_Celltype)[1], sort.by="none")$table
colnames(Result)<-c("Peak","logFC","logCPM","LR","PValue","FDR")
Result <- cbind(Result,TMM)

head(Result)
which(Result$Peak == "chr10_149311748_149312391")
print(sum(Result$FDR < 0.05))
#Result[Result$FDR < 0.05,][c(1:100),]
write.table(Result, file=paste0(CT,".EdgeRResult_PseudoReplicate_withPromoterRegion.txt"), quote=F, row.names=F, col.names=T, sep="\t")
            

##1043 --> IM-OC
## positive logFC means Bif3 is higher
## negative logFC means A619 is higher
