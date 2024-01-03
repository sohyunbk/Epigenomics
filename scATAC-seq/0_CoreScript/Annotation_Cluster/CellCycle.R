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
library(matrixStats)
library(Matrix)
library(gtools)
library(preprocessCore)
library(pheatmap) 
library(RColorBrewer)
######### Let's try heatmap #################
#meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.REF_CELLs.metadata.txt"

## Gene 
## Table shouldbe:
######## gene1 gene2 ... 
# cellA 1
# cellB 2

meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt" 
gene <- "/scratch/sb14489/3.scATAC/0.Data/CellCycle/CellCycle.bed"
GA <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619.txt"
AnnSlot <- "Ann_v4"
OutputPathFileName <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/2.CellCycle/A619_V4"


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
head(meta_data)
#### Calculate the average value in Cell type 
GroupInfo <- data.frame(barcode=rownames(meta_data),
                        Ann=meta_data[[AnnSlot]])
head(GroupInfo)

GBA_filtered <- merge(x = GroupInfo, y = GBA_filtered,
              by = "barcode", all.y = T)

dim(GBA_filtered)
head(gene_markers)

GroupInfo2<- data.frame(gene_name=gene_markers$geneID,
                        name=gene_markers$name)

head(GroupInfo2)
head(GBA_filtered)
GBA_filtered <- merge(x = GroupInfo2, y = GBA_filtered,
                      by = "gene_name", all.y = T)

head(GBA_filtered)
tail(GBA_filtered)
GBA_filtered$Ann_gene_name <- paste0(GBA_filtered$Ann,"&",GBA_filtered$gene_name)
head(GBA_filtered)

library(tidyverse)
MeanTable <- GBA_filtered %>%
  group_by(Ann_gene_name) %>%
  summarise_at(vars(accessability), list(MeanAcc = sum))
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

## Normalization --> zscore
#z <- t(as.matrix(scale(t(as.matrix(MeanTable_df)))))
MeanTable_df_T <- t(MeanTable_df)
head(MeanTable_df_T)
MeanTable_df_T <- cpm(MeanTable_df_T, log=TRUE, prior.count=5)
MeanTable_df_T_normalized <- normalize.quantiles(MeanTable_df_T)
Table_df_T_z <- t(as.matrix(scale(t(as.matrix(MeanTable_df_T_normalized)))))
head(Table_df_T_z)
colnames(Table_df_T_z) <- rownames(MeanTable_df)
rownames(Table_df_T_z) <- colnames(MeanTable_df)
z <- t(Table_df_T_z)
head(z)



#saveRDS(MeanTable_df, file="Heatmap_AverageAc.rds")
head(gene_markers)
anno <-data.frame(GeneID=gene_markers$geneID, group=gene_markers$name)
newCols <- colorRampPalette(grDevices::rainbow(length(unique(anno$group))))
head(anno)

annoCol <- newCols(length(unique(anno$group)))
head(z)

## Sort GeneID

head(gene_markers)
head(z)[,c(1:10)]
anno <- anno[order(anno$group),]
anno <- data.frame(row.names=anno$GeneID,group=anno$group)

z_sorted <- z[,rownames(anno)]
head(z_sorted)


#pdf(paste0(OutputPathFileName,"Heatmap_Average_z_CellCycle_All_group.pdf"), width =30, height=5, onefile=FALSE) 
#pheatmap(z_sorted,cluster_rows = F, annotation = anno,
#       cluster_cols = F) #, scale ="row" #cluster_cols = T # fontsize_row = 3.5,
#Heatmap(z_sorted,cluster_rows = F, 
#        top_annotation = HeatmapAnnotation(annotation = anno$group)) 
#dev.off() 

colnames(z_sorted)[c(1:10)]
head(anno)
#annotation <- data.frame(Var1 = factor(1:10 %% 2 == 0, labels = c("Exp1", "None")), Var1 = factor(1:10 %% 2 == 0, labels = c("Exp1", "None")))

#############################################
#############################################
###
## Draw BarPlot Ratio!
###
#############################################
#############################################

#colnames(MeanTable_df_T_normalized) <- rownames(MeanTable_df)
#rownames(MeanTable_df_T_normalized) <- colnames(MeanTable_df)
#head(MeanTable_df_T_normalized)
#head(anno)
BarPlotData <- data.frame()
for( i in levels(factor(anno$group))){
 GeneNames <- rownames(anno)[which(anno$group ==i)]
 Temp <- MeanTable_df_T[GeneNames,]
 Ratio <- colSums(Temp)/ colSums(MeanTable_df_T)
 value <- colSums(Temp)
 TempTable <- data.frame(value,Cellcycle=i)
 TempTable$Celltype <- rownames(TempTable)
 rownames(TempTable) <- NULL
 BarPlotData <- rbind(BarPlotData,TempTable)
 }
head(BarPlotData)

ggplot(BarPlotData, aes(fill=Cellcycle, y=value, x=Celltype)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(paste0(OutputPathFileName,"_stackedbarchart_CellCycle.pdf"), width=7, height=5)	


