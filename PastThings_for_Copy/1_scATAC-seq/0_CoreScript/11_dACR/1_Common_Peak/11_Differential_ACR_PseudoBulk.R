##### Back to EdgeR!! 
#conda activate r4-base
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)

## 1. Load files and get filtered Sparse file for save space
## --> should combine A619+Bif3 in the begining to keep all the peaks as features

Sparsefile_Bif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/Bif3_toComPeak.sparse"
Sparsefile_A619  <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/A619_toComPeak.sparse"
MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
InterGenic_peak <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/ComA619Bif3.unique500bpPeaks_BLRemove_Intergenic.bed")
dim(InterGenic_peak)

InterGenic_peak_Pos <- paste(InterGenic_peak$V1,InterGenic_peak$V2,InterGenic_peak$V3,sep="_")
head(InterGenic_peak_Pos)

GetFilteredSparseData <- function(Sparsefile,MetaFile,SelectedPeaksPos){
Sparse<- read.table(Sparsefile,header=F)
meta_data <- read.table(MetaFile, header=TRUE)
colnames(Sparse) <- c("PeakLocus","cellID","accessability")
## Remove!! chr6	181356755	181357256	1.10868807050461 !!!!!!!!!!!!!!
Sparse <- Sparse[Sparse$PeakLocus!="chr6_181356755_181357256",]

Sparse_Intergenic <- Sparse[Sparse$PeakLocus%in%SelectedPeaksPos,]
Sparse_Intergenic_SelectedCells <-Sparse_Intergenic[Sparse_Intergenic$cellID %in% meta_data$cellID,]
dim(Sparse_Intergenic_SelectedCells)
return(Sparse_Intergenic_SelectedCells)
}

Sparse_A619 <- GetFilteredSparseData(Sparsefile_A619,MetaFileA619,InterGenic_peak_Pos)
Sparse_Bif3 <- GetFilteredSparseData(Sparsefile_Bif3,MetaFileBif3,InterGenic_peak_Pos)
meta_data_A619 <- read.table(MetaFileA619, header=TRUE)
meta_data_Bif3 <- read.table(MetaFileBif3, header=TRUE)
CombinedMeta <- rbind(meta_data_A619,meta_data_Bif3)
head(CombinedMeta)

Sparse_Combined <- rbind(Sparse_A619,Sparse_Bif3)
Peak_Cell_Count <-  spread(Sparse_Combined,key = cellID,value =accessability)
dim(Peak_Cell_Count)
Peak_Cell_Count[is.na(Peak_Cell_Count)] <- 0

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_OnlyIntergenic")
saveRDS(Peak_Cell_Count, file="Peak_perCell_Counts.rds")
head(Peak_Cell_Count)[,c(1:10)]
#Peak_Cell_Count <- readRDS("Peak_perCell_Counts.rds")
Pos <- Peak_Cell_Count$PeakLocus


## list by cell types and Order the column of  count file 
cluster_name <- "Ann_v3"
celltypes <- levels(as.factor(CombinedMeta[,cluster_name]))
CountList <- list()
celltype <- celltypes[1]
for (celltype in celltypes){
  SelectedCellMeta <- CombinedMeta[which(CombinedMeta[,cluster_name] == celltype),]
  #head(SelectedCellMeta)
  SelectedCelltypeID <- rownames(SelectedCellMeta)
  Count_byCelltype <- Peak_Cell_Count[,SelectedCelltypeID]
  CountList[[celltype]] <- Count_byCelltype
}
saveRDS(CountList, file="CountMatrix_PerCelltype_list.rds")

#CountList <- readRDS("CountMatrix_PerCelltype_list.rds")
str(CountList)
### 2. EdgeR
celltypes <- names(CountList)

#### Make pseudobulk
Count_Celltype <- CountList[[celltypes[1]]]
dim(Count_Celltype)
head(Count_Celltype)[,c(1:3)]


DevideSampleIntoTwo <- function(StringVector){
  nSample <- length(StringVector)
  RandomSampleNumber <- sample(1:nSample, round(nSample/2), replace=F)
  Re1 <- StringVector[RandomSampleNumber]
  Re2 <- StringVector[!StringVector%in%Re1]
  return(list(Re1,Re2))
}
  
#https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# 3.3 Experiments with all combinations of multiple fac- tors

# Start with first cell type
ResultList <- list() ## But the loop might take too long so I will run it by submitting the jobs
nCT <-1
for (nCT in c(1:length(celltypes))){

CT<- celltypes[nCT]
print(CT)
Count_Celltype <- CountList[[CT]]
dim(Count_Celltype)
head(Count_Celltype)[,c(1:10)]
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

#TMM <- cpm(y, normalized.lib.sizes=TRUE,log=T)
#TMMcolName <- colnames(TMM)
#TMMcolName <- paste("TMM_",TMMcolName,sep="")
#str(TMMcolName)
#colnames(TMM) <- as.factor(TMMcolName)

lrt <- glmLRT(fit,coef=2)
Result <- topTags(lrt, n=dim(Count_Celltype)[1], sort.by="none")$table
colnames(Result)<-c("Peak","logFC","logCPM","LR","PValue","FDR")

head(Result)
print(sum(Result$FDR < 0.05))
#Result[Result$FDR < 0.05,][c(1:100),]
write.table(Result, file=paste0(CT,".EdgeRResult_PseudoReplicate.txt"), quote=F, row.names=F, col.names=T, sep="\t")
            
ResultList[[CT]] <- Result
}
##1043 --> IM-OC
LowFDR_IMOC <- ResultList[["IM-OC"]][ResultList[["IM-OC"]]$FDR < 0.05,]
head(LowFDR_IMOC)
sum(LowFDR_IMOC$logFC > 0) ## positive logFC means Bif3 is higher
sum(LowFDR_IMOC$logFC < 0) ## negative logFC means A619 is higher

saveRDS(ResultList, file="2Pesudo_Replicate_EdgeR_Result.rds")

sum(LowFDR_IMOC$logFC < -4)
LowFDR_IMOC[LowFDR_IMOC$logFC >3,]
LowFDR_IMOC[LowFDR_IMOC$logFC < -2,]

## Check distribution
length(Count_Celltype[1,])
ggplot() +aes(Count_Celltype[1,])+ geom_density()     
c(Count_Celltype[1,])
