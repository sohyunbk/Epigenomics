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
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,header=TRUE)

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
  return(nearest_genes_id)
}

DEGInfo_Bif3Higher <- DEGInfo[(DEGInfo$FDR < 0.05) & (DEGInfo$logFC > 0),]
DEGInfo_A619Higher <- DEGInfo[(DEGInfo$FDR < 0.05) & (DEGInfo$logFC < 0),]

SelectedPeak <- DEGInfo_Bif3Higher


meta <-"/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt"
gene <- "/scratch/sb14489/3.scATAC/0.Data/CellCycle/CellCycle.bed"
#GA <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_includingZmCLE7.txt"
#CellOrderFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt"
#OutputPathFileName <-"/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/2.CellCycle/A619"

meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)
GBA <- read.delim(GA)
CellOrder <- readLines(CellOrderFile)

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
#GroupInfo <- data.frame(barcode=rownames(meta_data),
#                        Ann=meta_data$Ann_V1)
GroupInfo <- data.frame(barcode=rownames(meta_data),
                        Ann=meta_data$Ann_v4)
head(GroupInfo)

GBA_filtered <- merge(x = GroupInfo, y = GBA_filtered,
                      by = "barcode", all.y = T)

dim(GBA_filtered)
head(gene_markers)

GroupInfo2<- data.frame(gene_name=gene_markers$geneID,
                        name=gene_markers$name)

head(GroupInfo2)
head(GBA_filtered)
dim(GBA_filtered)
GBA_filtered <- merge(x = GroupInfo2, y = GBA_filtered,
                      by = "gene_name", all.y = T)

head(GBA_filtered)
dim(GBA_filtered)

tail(GBA_filtered)
GBA_filtered$Ann_gene_name <- paste0(GBA_filtered$Ann,"&",GBA_filtered$gene_name)
head(GBA_filtered)

## It's meanTable but it's acutally sum table.
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
z <- as.matrix(scale(as.matrix(MeanTable_df)))
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

pdf(paste0(OutputPathFileName,"_Heatmap_Average_z_CellCycle_All_group.pdf"), width =30, height=5, onefile=FALSE) 
pheatmap(z_sorted,cluster_rows = F, annotation = anno,
         cluster_cols = F) #, scale ="row" #cluster_cols = T # fontsize_row = 3.5,
dev.off() 




########## By CellCycle ##

## Option 1: QunatileNormalization
QuantileNormalizedGenes_byCelltype <- normalize.quantiles(t(MeanTable_df))
MeanTable_df_T <- as.data.frame(QuantileNormalizedGenes_byCelltype)
colnames(MeanTable_df_T) <- rownames(MeanTable_df)
rownames(MeanTable_df_T) <- colnames(MeanTable_df)
head(MeanTable_df_T)

## Option 2: No QunatileNormalization
MeanTable_df_T <- t(MeanTable_df)
head(MeanTable_df_T)
#MeanTable_df_T <- cpm(MeanTable_df_T, log=TRUE, prior.count=5)
BarPlotData <- data.frame()

## ==> both option gives similar results

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


CellCycleColor <- c("#064f43","#4f4806","#4f0633","#38025e","#06274f")
BarPlotData$Celltype <- factor(BarPlotData$Celltype,levels=CellOrder)
BarPlotData$Cellcycle <- factor(BarPlotData$Cellcycle,
                                levels=c("G1","G1/S","S","G2/M","M"))
ggplot(BarPlotData, aes(fill=Cellcycle, y=value, x=Celltype)) + 
  scale_fill_manual(values = CellCycleColor) + # Replace with your colors
  geom_bar(position="fill", stat="identity")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(paste0(OutputPathFileName,"_stackedbarchart_CellCycle.pdf"), width=7, height=5)	


###################
## To write the statistic Table

BarPlotData_summary <- BarPlotData %>%
  group_by(Celltype, Cellcycle) %>%
  summarise(TotalValue = sum(value)) %>%
  ungroup()

# Calculate the total values per Celltype
TotalPerCelltype <- BarPlotData_summary %>%
  group_by(Celltype) %>%
  summarise(TotalValuePerCelltype = sum(TotalValue)) %>%
  ungroup()

# Join this with the original summary
BarPlotData_joined <- left_join(BarPlotData_summary, TotalPerCelltype, by = "Celltype")

# Calculate the ratio
BarPlotData_ratios <- BarPlotData_joined %>%
  mutate(Ratio = TotalValue / TotalValuePerCelltype) %>%
  select(Celltype, Cellcycle, Ratio)

# Spread the data so that each Cellcycle category becomes a column
BarPlotData_spread <- BarPlotData_ratios %>%
  spread(key = Cellcycle, value = Ratio, fill = 0)

# View the result
BarPlotData_spread
write.table(BarPlotData_spread, 
            paste0(OutputPathFileName,".CelltypeRatio.txt"),
            quote=F, row.names=F, col.names=T, sep="\t")
