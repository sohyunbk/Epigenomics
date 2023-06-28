## 202208 Got it from Pablo
## 202209 Edited by Sohyun
## Dotplot for all / part of markers
## conda activate r_env
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
######### Let's try heatmap #################
#meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.REF_CELLs.metadata.txt"

## Gene 
## Table shouldbe:
######## gene1 gene2 ... 
# cellA 1
# cellB 2

meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt" 
gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/221110_EarMarker.txt"
#GA <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/GA_A619.txt"
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)
GBA <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/MarkerGeneChange221110/opt_allgenes_impute.activity.rds")
head(GBA,n=10)
head(meta_data)
dim(gene_markers)
head(gene_markers)
dim(GBA)
GBA_filtered <- GBA[rownames(GBA) %in% gene_markers$geneID,]
dim(GBA_filtered)
dim(meta_data)
GBA_filtered <- GBA_filtered[,colnames(GBA_filtered) %in% rownames(meta_data)]
dim(GBA_filtered)
head(GBA_filtered)[,c(1:10)]

#### Calculate the average value in Cell type 
GroupInfo <- data.frame(barcode=rownames(meta_data),
                        Ann=meta_data$Ann_v3)
dim(GroupInfo)
head(GroupInfo)
GroupInfo <- GroupInfo[GroupInfo$barcode %in% colnames(GBA_filtered),]
Celltypes <- levels(factor(GroupInfo$Ann))

## Change gene name
head(GBA_filtered)[,c(1:10)]
rownames(gene_markers) <- gene_markers$geneID
dim(gene_markers)
rownames(GBA_filtered) <- gene_markers[rownames(GBA_filtered),]$name

SumTable <- matrix(0,nrow=nrow(gene_markers),ncol=length(Celltypes))

for(i in c(1:length(Celltypes))){
  SelectedB <- c(GroupInfo[GroupInfo$Ann==Celltypes[i],]$barcode)
  CellTypeValue <- GBA_filtered[,SelectedB]
  #SumTable[,i] <- Matrix::rowSums(CellTypeValue)
  SumTable[,i] <- Matrix::rowSums(CellTypeValue)
}
colnames(SumTable) <- as.factor(Celltypes)
rownames(SumTable) <- as.factor(rownames(GBA_filtered))
head(SumTable)
dim(SumTable)

library(matrixStats)
library(Matrix)
library(gtools)
library(preprocessCore)
#“average value in cluster * gene” --> CPM --> Zscore?
CPM <- cpm(SumTable, log=TRUE, prior.count=5)
#?normalize.quantiles
mat <- normalize.quantiles(CPM)
z <- t(as.matrix(scale(t(as.matrix(mat)))))
head(z)
colnames(z) <- as.factor(Celltypes)
rownames(z) <- as.factor(rownames(GBA_filtered))

#?cpm
#head(CPM)
#library(preprocessCore)
#?normalize.quantiles
#CPM <- normalize.quantiles(CPM)

Plotdata <- as.matrix(t(z))
head(Plotdata)

str(Plotdata)
typeof(Plotdata)
table(Plotdata)
library(pheatmap) 
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/2.HeatMap_Dotplot")

pdf("Heatmap_AnnotationV3_Imputation.pdf", width =15, height=5, onefile=FALSE) 
pheatmap(Plotdata, fontsize_row = 6,
         cluster_rows = TRUE) #, scale ="row" #cluster_cols = T
dev.off() 

library(ComplexHeatmap)
pdf("Heatmap_AnnotationV3_Imputation_complexheatmap.pdf", width =30, height=5, onefile=FALSE) 
Heatmap(matrix = as.matrix(Plotdata), 
        column_names_rot = 45,row_names_max_width = unit(15, "cm"),
        name = "norm")
dev.off() 

## Ex
library(ComplexHeatmap)
library(circlize)

set.seed(123)
mat = matrix(rnorm(80, 2), 8, 10)
mat = rbind(mat, matrix(rnorm(40, -2), 4, 10))
rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)

ha_column = HeatmapAnnotation(df = data.frame(type1 = c(rep("a", 5), rep("b", 5))),
                              col = list(type1 = c("a" =  "red", "b" = "blue")))
ha_row = rowAnnotation(df = data.frame(type2 = c(rep("A", 6), rep("B", 6))),
                       col = list(type2 = c("A" =  "green", "B" = "orange")), width = unit(1, "cm"))

ht1 = Heatmap(mat, name = "ht1", column_title = "Heatmap 1", top_annotation = ha_column)
ht2 = Heatmap(mat, name = "ht2", column_title = "Heatmap 2")
ht_list = ht1 + ht2 + ha_row

draw(ht_list)




