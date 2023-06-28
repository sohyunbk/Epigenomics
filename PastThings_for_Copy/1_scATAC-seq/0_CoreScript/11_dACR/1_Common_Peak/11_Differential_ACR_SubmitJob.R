##### Back to EdgeR!! 
#conda activate r4-base
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)

args <- commandArgs(T)
CT <- as.character(args[1])

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_OnlyIntergenic")

Peak_Cell_Count <- readRDS("Peak_perCell_Counts.rds")
Pos <- Peak_Cell_Count$PeakLocus
CountList <- readRDS("CountMatrix_PerCelltype_list.rds")
celltypes <- names(CountList)

# Start with first cell type
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
#cbind(targets,Group=Group)
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
colnames(Result)<-c("Peak","logFC","logCPM","LR","PValue","FDR")

head(Result)
print(sum(Result$FDR < 0.05))
#Result[Result$FDR < 0.05,][c(1:100),]
write.table(Result, file=paste0(CT,".EdgeRResult.txt"), quote=F, row.names=F, col.names=T, sep="\t")
            
##1043 --> IM-OC

