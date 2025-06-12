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
## load gene*cell table.
WT_GeneXCT <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC/A619_AnnV4.GeneBodyACC.byGeneXCT.txt",header=TRUE)
Bif3_GeneXCT <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC/Bif3_AnnV4.GeneBodyACC.byGeneXCT.txt",header=TRUE)
colnames(Bif3_GeneXCT)[-1] <- paste("Bif3&", colnames(Bif3_GeneXCT)[-1], sep="")
colnames(WT_GeneXCT)[-1] <- paste("A619&", colnames(WT_GeneXCT)[-1], sep="")
head(Bif3_GeneXCT)
head(WT_GeneXCT)
GeneXCT <- merge(WT_GeneXCT, Bif3_GeneXCT, by = "gene")
head(GeneXCT)
## normalization
GeneNames <- GeneXCT[, 1]
gene_counts_df <- GeneXCT[, -1]

## To filter some genes with very low Tn5.
library(ggplot2)
library(tidyr)
# Reshape the data into long format
gene_counts_long <- gather(gene_counts_df, key = "gene", value = "count")
# Create a density plot
ggplot(gene_counts_long, aes(x = count)) +
  geom_density(fill = "skyblue", color = "blue") +
  labs(title = "Density Plot of Gene Counts", x = "Count", y = "Density")+
  scale_x_continuous(breaks = seq(0, 1000, by = 50))
ggsave("DensityPlot_Tn5_GeneBodyACC.pdf" , width=40, height=5)
## Set up the cut off to 50!
head(GeneXCT)
library(edgeR)
library(preprocessCore)
GeneXCT_MoreThan50Tn5 <- GeneXCT %>%
  filter_all(all_vars(. > 50))
dim(GeneXCT_MoreThan50Tn5)
GeneNames <- GeneXCT_MoreThan50Tn5[, 1]
gene_counts_df <- GeneXCT_MoreThan50Tn5[, -1]

cpm_data <- cpm(DGEList(counts = gene_counts_df), log = FALSE)
qnorm_data <- normalize.quantiles(as.matrix(gene_counts_df))
rownames(qnorm_data) <- GeneNames
colnames(qnorm_data) <- colnames(GeneXCT[,-1])
head(qnorm_data)
dim(qnorm_data)
#write.table(qnorm_data, file=paste0(SampleName,".GeneBodyACC.byGeneXCT.CPMQuantilNor.txt"), quote=F, row.names=T, col.names=T, sep="\t")
##############
head(Bif3HigherTAAT_Nearest_SecondN)
head(Bif3HigherTAAT_Nearest_SecondN$Nearest_Gene)
length(Bif3HigherTAAT_Nearest_SecondN$Nearest_Gene)

CTOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt")
CTOrder <- gsub("-", ".", CTOrder)
new_CTOrder <- unlist(lapply(CTOrder, function(x) c(paste("A619&", x, sep=""), paste("Bif3&", x, sep=""))))
head(new_CTOrder)

dim(qnorm_data)
qnorm_data_orderd <- qnorm_data[,new_CTOrder]
QNorm_SelectedGene <- qnorm_data_orderd[rownames(qnorm_data_orderd)%in%Bif3HigherTAAT_Nearest_SecondN$Nearest_Gene,]
dim(QNorm_SelectedGene)
head(QNorm_SelectedGene)
library(gplots)

# Define color palette for the heatmap
my_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(255)
# Create heatmap
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC")


###### logFC of gene bodyACC
FCTable <- data.frame(
  L1 = QNorm_SelectedGene[,"Bif3&L1"] / QNorm_SelectedGene[,"A619&L1"],
  L1atFloralMeristem = QNorm_SelectedGene[,"Bif3&L1atFloralMeristem"] / QNorm_SelectedGene[,"A619&L1atFloralMeristem"],
  FloralMeristem_SuppressedBract = QNorm_SelectedGene[,"Bif3&FloralMeristem_SuppressedBract"] / QNorm_SelectedGene[,"A619&FloralMeristem_SuppressedBract"],
  IM_OC = QNorm_SelectedGene[,"Bif3&IM.OC"] / QNorm_SelectedGene[,"A619&IM.OC"],
  SPM_base_SM_base = QNorm_SelectedGene[,"Bif3&SPM.base_SM.base"] / QNorm_SelectedGene[,"A619&SPM.base_SM.base"],
  ProcambialMeristem_ProtoXylem_MetaXylem = QNorm_SelectedGene[,"Bif3&ProcambialMeristem_ProtoXylem_MetaXylem"] / QNorm_SelectedGene[,"A619&ProcambialMeristem_ProtoXylem_MetaXylem"],
  PhloemPrecursor = QNorm_SelectedGene[,"Bif3&PhloemPrecursor"] / QNorm_SelectedGene[,"A619&PhloemPrecursor"],
  ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma = QNorm_SelectedGene[,"Bif3&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma"] / QNorm_SelectedGene[,"A619&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma"],
  XylemParenchyma_PithParenchyma = QNorm_SelectedGene[,"Bif3&XylemParenchyma_PithParenchyma"] / QNorm_SelectedGene[,"A619&XylemParenchyma_PithParenchyma"],
  Unknown1 = QNorm_SelectedGene[,"Bif3&Unknown1"] / QNorm_SelectedGene[,"A619&Unknown1"],
  Unknown2 = QNorm_SelectedGene[,"Bif3&Unknown2"] / QNorm_SelectedGene[,"A619&Unknown2"],
  Unknown_Sclerenchyma = QNorm_SelectedGene[,"Bif3&Unknown_Sclerenchyma"] / QNorm_SelectedGene[,"A619&Unknown_Sclerenchyma"],
  Unknown_lowFRiP = QNorm_SelectedGene[,"Bif3&Unknown_lowFRiP"] / QNorm_SelectedGene[,"A619&Unknown_lowFRiP"],
  G2_M = QNorm_SelectedGene[,"Bif3&G2_M"] / QNorm_SelectedGene[,"A619&G2_M"]
)
head(FCTable)
dim(FCTable)
FCTable_matrix <- as.matrix(log2(FCTable))
head(FCTable_matrix)
IM_OC_column <- FCTable_matrix[,"IM_OC"]
ordered_rows <- order(IM_OC_column)
FCTable_ordered <- FCTable_matrix[ordered_rows, ]
head(FCTable_ordered)
tail(FCTable_ordered)

### Let's match it with the gene symbol
GeneInfo <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.Curated.txt", fill = TRUE,header=TRUE)
head(GeneInfo)
GeneID_GeneSymbol <- data.frame(GeneInfo$gene_model,GeneInfo$locus_symbol)
head(GeneID_GeneSymbol)
head(FCTable_ordered)
FCTable_ordered_geneSymbol <- merge(GeneID_GeneSymbol, FCTable_ordered, by.x = "GeneInfo.gene_model", by.y = "row.names", all.y = TRUE)
head(FCTable_ordered_geneSymbol)

FCTable_ordered_geneSymbol$GeneInfo.gene_model <- ifelse(is.na(FCTable_ordered_geneSymbol$GeneInfo.locus_symbol),
                                                         FCTable_ordered_geneSymbol$GeneInfo.gene_model,
                                                         FCTable_ordered_geneSymbol$GeneInfo.locus_symbol)
head(FCTable_ordered_geneSymbol)
# Remove GeneInfo.locus_symbol column
FCTable_ordered_geneSymbol$GeneInfo.locus_symbol <- NULL

head(FCTable_ordered_geneSymbol)
rownames(FCTable_ordered_geneSymbol) <- FCTable_ordered_geneSymbol[,1]
FCTable_ordered_geneSymbol <- FCTable_ordered_geneSymbol[,-1]
head(FCTable_ordered_geneSymbol)
ordered_rows <- order(FCTable_ordered_geneSymbol[,"IM_OC"])
FCTable_ordered_geneSymbol <- FCTable_ordered_geneSymbol[ordered_rows, ]
head(FCTable_ordered_geneSymbol)
FCTable_ordered_geneSymbol <- as.matrix(FCTable_ordered_geneSymbol)

my_palette <- colorRampPalette(c("blue", "white", "red"))(101)  # Adjusted to 101 colors
breaks <- c(-2.1, seq(-2, 2, length.out = 100), 2.1)  # Adjusted breaks
library(fields)

pdf("logFC_TAATMotifdACRClosetGene_HeatMap.pdf", width = 10, height = 8)  # Specify the file name and dimensions

# Create the heatmap
heatmap(FCTable_ordered_geneSymbol,
        col = my_palette,
        breaks = breaks,
        scale = "none",
        margins = c(20,20),  # Increase the left margin to accommodate longer column labels
        main = "Heatmap of QNorm_SelectedGene",
        Rowv = NA,
        Colv = NA,
        cexCol = 1)  # Adjust the size of the column labels if necessary

# Adjust margins to make space for the legend
par(mar = c(5, 4, 2, 6))  # Reduce right margin to 6

# Add color legend manually with custom size and position, and specific values
image.plot(zlim = range(breaks),
           col = my_palette,
           legend.only = TRUE,
           horizontal = FALSE,
           axis.args = list(at = c(-2.1, 0, 2.1), labels = c(-2.1, 0, 2.1), cex.axis = 0.5),  # Specify values for the legend
           legend.width = 0.8,  # Adjust the width of the legend box
           legend.shrink = 0.5
           )  # Adjust the size of the color legend bar

dev.off()

pdf("logFC_TAATMotifdACRClosetGene_HeatMap_rowAutoOrder.pdf", width = 10, height = 8)  # Specify the file name and dimensions

# Create the heatmap
heatmap(FCTable_ordered_geneSymbol,
        col = my_palette,
        breaks = breaks,
        scale = "none",
        margins = c(20,20),  # Increase the left margin to accommodate longer column labels
        main = "Heatmap of QNorm_SelectedGene",
        Rowv = TRUE,
        Colv = NA,
        cexCol = 1)  # Adjust the size of the column labels if necessary

# Adjust margins to make space for the legend
par(mar = c(5, 4, 2, 6))  # Reduce right margin to 6

# Add color legend manually with custom size and position, and specific values
image.plot(zlim = range(breaks),
           col = my_palette,
           legend.only = TRUE,
           horizontal = FALSE,
           axis.args = list(at = c(-2.1, 0, 2.1), labels = c(-2.1, 0, 2.1), cex.axis = 0.5),  # Specify values for the legend
           legend.width = 0.8,  # Adjust the width of the legend box
           legend.shrink = 0.5
)  # Adjust the size of the color legend bar

dev.off()

write.table(FCTable_ordered_geneSymbol,"logFCNeartGene_TAATdACR.txt", 
            quote=F, row.names=F, col.names=T, sep="\t")


############################################
############# ARF gene!! ###################
## 1) Get all the ARF gene.