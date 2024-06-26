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
######### Let's try heatmap #################
#meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.REF_CELLs.metadata.txt"

## Gene 
## Table shouldbe:
######## gene1 gene2 ... 
# cellA 1
# cellB 2

meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_withSub_withName.metadata.txt" 
gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/Ear.markers.v5.InSitu.All.Filtered.txt"
GA <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/GA_A619.txt"

meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)
GBA <- read.delim(GA)

head(GBA)
head(meta_data)
dim(gene_markers)
head(gene_markers)
dim(GBA)
GBA_filtered <- GBA[GBA$gene_name %in% gene_markers$geneID,]
dim(GBA_filtered)
GBA_filtered <- GBA_filtered[GBA_filtered$barcode %in% rownames(meta_data),]
dim(GBA_filtered)
head(GBA_filtered)

#### Calculate the average value in Cell type 
GroupInfo <- data.frame(barcode=rownames(meta_data),
                        Ann=meta_data$Ann_V1)
head(GroupInfo)
GBA_filtered <- merge(x = GroupInfo, y = GBA_filtered,
              by = "barcode", all.y = T)
dim(GBA_filtered)
GroupInfo2<- data.frame(gene_name=gene_markers$geneID,
                        name=gene_markers$name)
GBA_filtered <- merge(x = GroupInfo2, y = GBA_filtered,
                      by = "gene_name", all.y = T)

head(GBA_filtered)
tail(GBA_filtered)
GBA_filtered$Ann_gene_name <- paste0(GBA_filtered$Ann,"&",GBA_filtered$name)
head(GBA_filtered)

library(tidyverse)
MeanTable <- GBA_filtered %>%
  group_by(Ann_gene_name) %>%
  summarise_at(vars(accessability), list(MeanAcc = mean))
head(MeanTable)
MeanTable<- data.frame(celltype = as.character(lapply(strsplit(as.character(MeanTable$Ann_gene_name),
                                                                split="&"), "[", 1)),
                       gene= as.character(lapply(strsplit(as.character(MeanTable$Ann_gene_name),
                                                             split="&"), "[", 2)),
                       accessability = as.character(MeanTable$MeanAcc))
head(MeanTable)
levels(factor(MeanTable$celltype))
str(MeanTable)
MeanTable$accessability <- as.numeric(MeanTable$accessability) 

## Spread a pair of columns into a field of cells
library(tidyverse)
#https://rstudio-education.github.io/tidyverse-cookbook/tidy.html
MeanTable_spread <- as_tibble(MeanTable)
head(MeanTable_spread)
 
MeanTable_spread<- MeanTable_spread %>% 
  spread(key = gene, value = accessability,fill=0)
head(MeanTable_spread)
dim(MeanTable_spread)


#GBA_filtered 
#GBA_filtered[is.na(GBA_filtered)] <- 0

MeanTable_df <- as.data.frame(MeanTable_spread)
head(MeanTable_df)
rownames(MeanTable_df) <- MeanTable_df$celltype

head(MeanTable_df)[,c(1:10)]
MeanTable_df <- subset(MeanTable_df, select=-c(celltype))
tail(MeanTable_df)
dim(MeanTable_df)
MeanTable_df[is.na(MeanTable_df)]
str(MeanTable_df)
#as.numeric(MeanTable_df)

library(pheatmap) 
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/2.HeatMap_Dotplot")

saveRDS(MeanTable_df, file="Heatmap_AverageAc.rds")

pdf("Heatmap_AverageAc.pdf", width =15, height=5, onefile=FALSE) 
pheatmap(MeanTable_df,cluster_rows = F , fontsize_row = 3.5,
         scale="column",cluster_cols = T) #, scale ="row" #cluster_cols = T
dev.off() 



