## 202208 Got it from Pablo
## 202209 Edited by Sohyun
## Dotplot for all / part of markers

# load libraries
library(dplyr)
library("edgeR")
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(tidyverse)
library(tidyr)
library(pheatmap) 
library(RColorBrewer)
library("optparse")
library(preprocessCore)
library(devtools)
library("fgsea")
library("here")
library(devtools)
library(tidyverse)
library(Matrix)
library(magrittr) # needs to be run every time you start R and want to use %>%
library("optparse")
library(GenomicRanges)
library(ggplot2)

#### 1) Find the two closest genes!
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/dACR_withTAATInfo.txt"
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,fill=TRUE,header=TRUE)

## Function1: Get nearest genes by cutoff of dACR.
Get_NearestGenes <-function(SelectedPeak){
  genes_ranges <- GRanges(seqnames=genes_data$V1, 
                          ranges=IRanges(start=genes_data$V2, end=genes_data$V3),
                          gene_id=genes_data$V4)
  #head(SelectedPeak)
  split_names <- strsplit(SelectedPeak$Peak, split = "_")
  chromosomes <- sapply(split_names, function(x) x[1])
  starts <- as.integer(sapply(split_names, function(x) x[2]))
  ends <- as.integer(sapply(split_names, function(x) x[3]))
  granges_object <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))
  results <- data.frame(ACR = character(), Nearest_Gene = character(), Second_Nearest_Gene = character(), stringsAsFactors = FALSE)
  for (i in seq_along(granges_object)) {
    single_acr <- granges_object[i]
    # Find the nearest gene
    nearest_gene <- distanceToNearest(single_acr, genes_ranges)
    nearest_index <- subjectHits(nearest_gene)
    # Exclude the nearest gene and find the second nearest
    genes_minus_closest <- genes_ranges[-nearest_index]
    second_nearest_gene <- distanceToNearest(single_acr, genes_minus_closest)
    second_nearest_index <- subjectHits(second_nearest_gene)
    # Store the results
    results <- rbind(results, data.frame(
      ACR = as.character(SelectedPeak$Peak[i]),
      Nearest_Gene = as.character(mcols(genes_ranges[nearest_index])$gene_id),
      Second_Nearest_Gene = as.character(mcols(genes_minus_closest[second_nearest_index])$gene_id)
    ))
  }
  return(results)
}

DEGInfo_Bif3Higher <- DEGInfo[(DEGInfo$FDR < 0.05) & (DEGInfo$logFC > 0),]
DEGInfo_Bif3Higher_TAAT <- DEGInfo[(DEGInfo$FDR < 0.05) & (DEGInfo$logFC > 0) &
                                     (DEGInfo$TAAT == "TAAT"),]

DEGInfo_A619Higher <- DEGInfo[(DEGInfo$FDR < 0.05) & (DEGInfo$logFC < 0),]
dim(DEGInfo_Bif3Higher_TAAT)

Bif3HigherTAAT_Nearest_SecondN <- Get_NearestGenes(DEGInfo_Bif3Higher_TAAT)
NearestGenes <- c(Bif3HigherTAAT_Nearest_SecondN$Nearest_Gene,
                  Bif3HigherTAAT_Nearest_SecondN$Second_Nearest_Gene)


### 2) Get Gene body acc for all the genes 


GBA <- read.delim(GA)
meta_data <- read.delim(meta)
GBA_filtered <- GBA[GBA$barcode %in% rownames(meta_data),]
GroupInfo <- data.frame(barcode=rownames(meta_data),
                        Ann=meta_data$Ann_v4)
head(GroupInfo)

GBA_filtered <- merge(x = GroupInfo, y = GBA_filtered,
                      by = "barcode", all.y = T)

dim(GBA_filtered)
head(gene_markers)
head(GBA_filtered)

GBA_filtered$Ann_gene_name <- paste0(GBA_filtered$Ann,"&",GBA_filtered$gene_name)
head(GBA_filtered)

SumTable <- GBA_filtered %>%
  group_by(Ann_gene_name) %>%
  summarise_at(vars(accessability), list(SumAcc = sum))
head(SumTable)

SumTable<- data.frame(celltype = as.character(lapply(strsplit(as.character(SumTable$Ann_gene_name),
                                                               split="&"), "[", 1)),
                       gene= as.character(lapply(strsplit(as.character(SumTable$Ann_gene_name),
                                                          split="&"), "[", 2)),
                       accessability = as.character(SumTable$SumAcc))
head(SumTable)
levels(factor(SumTable$celltype))
str(SumTable)
SumTable$accessability <- as.numeric(SumTable$accessability) 
SumTable_spread <- as_tibble(SumTable)
head(SumTable_spread)
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC/")
SampleName <-"A619"
SumTable_spread<- SumTable_spread %>% 
  spread(key = celltype, value = accessability,fill=0)
head(SumTable_spread)
dim(SumTable_spread)
write.table(SumTable_spread, file=paste0(SampleName,".GeneBodyACC.byGeneXCT.txt"), quote=F, row.names=F, col.names=T, sep="\t")


## normalization
gene_counts_df <- as.data.frame(SumTable_spread)
rownames(gene_counts_df) <- gene_counts_df[, 1]
gene_counts_df <- gene_counts_df[, -1]
library(edgeR)
library(preprocessCore)

cpm_data <- cpm(DGEList(counts = gene_counts_df), log = FALSE)
qnorm_data <- normalize.quantiles(as.matrix(gene_counts_df))
rownames(qnorm_data) <- gene_counts_df[, 1]
colnames(qnorm_data) <- colnames(gene_counts_df)
head(qnorm_data)
write.table(qnorm_data, file=paste0(SampleName,".GeneBodyACC.byGeneXCT.CPMQuantilNor.txt"), quote=F, row.names=T, col.names=T, sep="\t")
##############


