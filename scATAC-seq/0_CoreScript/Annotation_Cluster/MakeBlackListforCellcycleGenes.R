library(dplyr)
library(here)
library(DESeq2)
library(tidyverse)
library(rlang)
library(GenomicRanges)
library("here")
library(devtools)
library(Seurat)
library(cicero)
library(Gviz)
library(rlang)
library("here")
library(devtools)
library(Seurat)
library(Matrix)
library(irlba)
#install.packages("Matrix", INSTALL_opts = '--no-lock')
#install.packages("latticeExtra", dependencies = TRUE, INSTALL_opts = '--no-lock')
load_all('/home/sb14489/Socrates')

## Make blacklsit again for cellcycle genes
obj <- list()

ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf"
attribute="Parent"
obj$gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format="gtf", dbxrefTag=attribute)))
obj$gff
  

FeatureName <- "gene"
## Select Features
if (FeatureName == "gene"){
  Ann <- genes(obj$gff)
} else if(FeatureName == "transcript"){
  Ann <- transcripts(obj$gff)
} else{message("ERROR: Feature name should be 'gene' or 'transcript")}


CellCycle <- "/scratch/sb14489/3.scATAC/0.Data/CellCycle/B73_v5_cell_cycle_genes.txt"
CellCycle <- read.table(CellCycle, header=FALSE)
head(CellCycle)
CellCycleGenes <- CellCycle$V1
head(BroadRange_Ann)

gr <- BroadRange_Ann[names(BroadRange_Ann)%in%CellCycleGenes,]
cc_bed <- data.frame(chr=seqnames(gr),
                                   start=start(gr),
                                   end=end(gr),
                                   geneID=names(gr))
head(cc_bed)
cc_bed$geneID 
CellCycle_Bed <- merge(cc_bed,CellCycle, by.x="geneID", by.y="V1")
CellCycle_Bed$type <- CellCycle_Bed$V2
col(CellCycle_Bed)
CellCycle_Bed <- CellCycle_Bed[,c(2,3,4,1,5,6)]
colnames(CellCycle_Bed) <- c("chr","start","end","geneID","name","type")
tail(CellCycle_Bed)
setwd("/scratch/sb14489/3.scATAC/0.Data/CellCycle")
write.table(CellCycle_Bed, file="CellCycle_UpDownStream500bp_withHeader.bed",
           quote=F, row.names=F, col.names=T, sep="\t")
#data.frame(v1=d$v1, v4=m[match(d$v2, m$v3), 2])
