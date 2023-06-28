##### Back to EdgeR!! 
#conda activate r4-base
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)

## 1. Load files and get filtered Sparse file for save space
## --> should combine A619+Bif3 in the begining to keep all the peaks as features

Sparsefile_Bif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/bif3_CommonpeakwithA619.sparse"
Sparsefile_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/A619_CommonpeakwithBif3.sparse"
MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
InterGenic_peak <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/ComA619Bif3.unique500bpPeaks_BLRemove_Intergenic.bed",header=F)
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

CountList <- readRDS("CountMatrix_PerCelltype_list.rds")
str(CountList)
### 2. EdgeR
celltypes <- names(CountList)

#https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# 3.3 Experiments with all combinations of multiple fac- tors

# Start with first cell type
ResultList <- list() ## But the loop might take too long so I will run it by submitting the jobs
nCT <-1
for (nCT in c(2:length(celltypes))){

CT<- celltypes[nCT]
print(CT)
Count_Celltype <- CountList[[CT]]
dim(Count_Celltype)
head(Count_Celltype)[,c(1:10)]

CellNames <- colnames(Count_Celltype)
y <- DGEList(counts=Count_Celltype, gene=Pos)
y <- calcNormFactors(y)

SampleName <- c(rep("A619",length(CellNames[endsWith(CellNames, 'A619')])))
SampleName <- c(SampleName,rep("A619",length(CellNames[endsWith(CellNames, 'A619_2')])))
SampleName <- c(SampleName,rep("bif3",length(CellNames[endsWith(CellNames, 'bif3')])))
SampleName <- c(SampleName,rep("bif3",length(CellNames[endsWith(CellNames, 'bif3_2')])))

Library <- c(rep("A619",length(CellNames[endsWith(CellNames, 'A619')])))
Library <- c(Library,rep("A619_2",length(CellNames[endsWith(CellNames, 'A619_2')])))
Library <- c(Library,rep("bif3",length(CellNames[endsWith(CellNames, 'bif3')])))
Library <- c(Library,rep("bif3_2",length(CellNames[endsWith(CellNames, 'bif3_2')])))

targets<-  data.frame(SampleName,Library)
Group <- factor(paste(SampleName,Library,sep="."))
cbind(targets,Group=Group)
#design <- model.matrix(~0+Group)
#colnames(design) <- levels(Group)
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
colnames(Result)<-c("Taxa","logFC","logCPM","LR","PValue","FDR")

head(Result)
print(sum(Result$FDR < 0.05))
#Result[Result$FDR < 0.05,][c(1:100),]
write.table(Result, file=paste0(CT,".EdgeRResult.txt"), quote=F, row.names=F, col.names=T, sep="\t")
            
ResultList[[CT]] <- Result
}
##1043 --> IM-OC

## Check distribution
length(Count_Celltype[1,])
ggplot() +aes(Count_Celltype[1,])+ geom_density()     
c(Count_Celltype[1,])
