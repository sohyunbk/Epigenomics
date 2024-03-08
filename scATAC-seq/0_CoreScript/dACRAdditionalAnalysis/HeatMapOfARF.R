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


### 1) Get Gene body acc for all the genes 
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

## Set up the cut off to 50!
head(GeneXCT)
dim(GeneXCT)
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
#### 2) Find ARF genes
GeneInfo <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.Curated.txt",
                       fill = TRUE,header=TRUE)
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/dACR_withTAATInfo.txt"
DEGInfo <- read.table(DEGFile,fill=TRUE,header=TRUE)
#head(GeneInfo)
GeneID_GeneSymbol <- data.frame(GeneInfo$gene_model,GeneInfo$locus_symbol)
head(GeneID_GeneSymbol)
GeneID_GeneSymbol_arf <- GeneID_GeneSymbol[grep("^arf", GeneID_GeneSymbol$GeneInfo.locus_symbol, 
                                                ignore.case = TRUE), ]
dim(GeneID_GeneSymbol_arf)

GeneID_GeneSymbol_arf <- GeneID_GeneSymbol[grep("^arf", GeneID_GeneSymbol$GeneInfo.locus_symbol, 
                                                ignore.case = TRUE), ]

DEGInfo_arf <- DEGInfo[grep("^arf", DEGInfo$locus_symbol, 
                                      ignore.case = TRUE), ]
dim(DEGInfo_arf)
head(DEGInfo_arf)
DEGInfo_arf$logFC
DEGInfo_arf$FDR
## 14 ACRs near ARFs
ARF_ACR <- data.frame(ACR = DEGInfo_arf$Peak,GeneName =DEGInfo_arf$locus_symbol,
                      geneid =DEGInfo_arf$gene_model,
                      ACR_logFC = DEGInfo_arf$logFC)
head()

head(qnorm_data)


SelectedGenes <- ARF_ACR$geneid
CTOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt")
CTOrder <- gsub("-", ".", CTOrder)
new_CTOrder <- unlist(lapply(CTOrder, function(x) c(paste("A619&", x, sep=""), paste("Bif3&", x, sep=""))))
head(new_CTOrder)

dim(qnorm_data)
qnorm_data_orderd <- qnorm_data[,new_CTOrder]

QNorm_SelectedGene <- qnorm_data_orderd[rownames(qnorm_data_orderd)%in%SelectedGenes,]
dim(QNorm_SelectedGene)
head(QNorm_SelectedGene)

library(gplots)

# Define color palette for the heatmap
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
min(FCTable_ordered)
dim(FCTable_ordered)
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
FCTable_ordered_geneSymbol
ARF_ACR
ARF_ACR_Combined <- data.frame(ARFName=NA,
                               ACR1_ACR=NA,
                               ACR1_ACR_logFC=NA,
                               ACR2_ACR=NA,
                               ACR2_ACR_logFC=NA
                          )
for (ARFName in unique(ARF_ACR$GeneName)){
  print(ARFName)
  Rows <- ARF_ACR[ARF_ACR$GeneName==ARFName,]
  if (nrow(Rows) ==2){
    Table1 <- Rows[1,-c(2,3)]
    colnames(Table1) <- paste0("ACR1","_",colnames(Table1)) 
    Table2 <- Rows[2,-c(2,3)]
    colnames(Table2) <- paste0("ACR2","_",colnames(Table2)) 
    NewRow <- cbind(ARFName = ARFName,Table1,Table2)
    ARF_ACR_Combined <- rbind(ARF_ACR_Combined,NewRow)
  }
  if (nrow(Rows)==1){
    Table1 <- Rows[1,-c(2,3)]
    colnames(Table1) <- paste0("ACR1","_",colnames(Table1)) 
    Table2 <- data.frame(ACR2_ACR=NA,
                         ACR2_ACR_logFC=0)
    
    NewRow <- cbind(ARFName = ARFName,Table1,Table2)
    ARF_ACR_Combined <- rbind(ARF_ACR_Combined,NewRow)
  }
}

ARF_ACR_Combined<- ARF_ACR_Combined[-1,]
head(ARF_ACR_Combined)
FCTable_ordered_geneSymbol
FCTable_ordered_geneSymbol$ARFName <- rownames(FCTable_ordered_geneSymbol)
ARF_all <- merge(ARF_ACR_Combined, FCTable_ordered_geneSymbol, by = "ARFName", all.x = TRUE)
rownames(ARF_all) <- ARF_all$ARFName
ARF_all <- ARF_all[rev(c("arftf4","arftf25","arftf30","arftf18","arftf23",
                     "arftf20","arftf26","arftf36","arftf3","arftf10")),]
write.table(ARF_all,"ARF_logFC.txt",
            quote=F, row.names=F, col.names=T, sep="\t")
ARF_GeneBody <- ARF_all[, !names(ARF_all) %in% c("ARFName", "ACR1_ACR", "ACR2_ACR","ACR1_ACR_logFC", "ACR2_ACR_logFC")]
ARF_dACR <- ARF_all[,  c("ACR1_ACR_logFC", "ACR2_ACR_logFC")]
ARF_dACR_rounded <- round(ARF_dACR, digits = 3)
ARF_dACR_rounded_matrix <- as.matrix(ARF_dACR_rounded)
ARF_GeneBody<-as.matrix(ARF_GeneBody)

breaks <- seq(-max(abs(ARF_dACR_rounded_matrix)), 
              max(abs(ARF_dACR_rounded_matrix)), length.out = 103)
my_palette <- colorRampPalette(c("#560e8f", "white", "red"))(102)  # Adjusted to 102 colors

pdf("ARF_dACR_Heatmap.pdf", width = 8, height = 12)  # Specify the file name and dimensions
# Create the heatmap with annotation
heatmap(ARF_dACR_rounded_matrix,
        col = my_palette,
        breaks = breaks, # Specify color palette
        scale = "none",
        main = "Heatmap of ARF_dACR",
        Rowv = NA, Colv = NA,
        cexRow = 0.5, cexCol = 0.8,  # Adjust row and column label sizes
        margins = c(10, 10))  # Adjust margins

dev.off()
library(fields)

pdf("ARF_dACR_ColorLegend.pdf", width = 4, height = 6)  # Specify the file name and dimensions
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", asp = 1, xlab = "", ylab = "")

# Add the color legend
image.plot(zlim = range(breaks),
           col = my_palette,
           legend.only = TRUE,
           horizontal = FALSE,
           axis.args = list(at = c(-max(abs(ARF_dACR_rounded_matrix)), 0, max(abs(ARF_dACR_rounded_matrix))), 
                            labels = c(-max(abs(ARF_dACR_rounded_matrix)), 0, max(abs(ARF_dACR_rounded_matrix))),
                            cex.axis = 0.5),  # Specify values for the legend
           legend.width = 0.8,  # Adjust the width of the legend box
           legend.shrink = 0.5)  # Adjust the size of the color legend bar

# Save the color legend plot to PDF
dev.off()

###
cols_to_keep <- !grepl("Unknown|G2", colnames(ARF_GeneBody))
ARF_GeneBody_NotUnknown <- ARF_GeneBody[, cols_to_keep]

breaks <- seq(-max(abs(ARF_GeneBody_NotUnknown)), 
              max(abs(ARF_GeneBody_NotUnknown)), length.out = 103)
my_palette <- colorRampPalette(c("blue", "white", "red"))(102)  # Adjusted to 102 colors



pdf("ARF_GenebodyACC_Heatmap.pdf", width = 30, height = 12)  # Specify the file name and dimensions
# Create the heatmap with annotation
heatmap(ARF_GeneBody_NotUnknown,
        col = my_palette,
        breaks = breaks, # Specify color palette
        scale = "none",
        main = "Heatmap of ARF_dACR",
        Rowv = NA, Colv = NA,
        cexRow = 0.5, cexCol = 0.8,  # Adjust row and column label sizes
        margins = c(10, 10))  # Adjust margins

dev.off()


pdf("ARF_GeneACC_ColorLegend.pdf", width = 4, height = 6)  # Specify the file name and dimensions
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i", asp = 1, xlab = "", ylab = "")

# Add the color legend
image.plot(zlim = range(breaks),
           col = my_palette,
           legend.only = TRUE,
           horizontal = FALSE,
           axis.args = list(at = c(-max(abs(ARF_GeneBody_NotUnknown)), 0, max(abs(ARF_GeneBody_NotUnknown))), 
                            labels = c(-max(abs(ARF_GeneBody_NotUnknown)), 0, max(abs(ARF_GeneBody_NotUnknown))),
                            cex.axis = 0.5),  # Specify values for the legend
           legend.width = 0.8,  # Adjust the width of the legend box
           legend.shrink = 0.5)  # Adjust the size of the color legend bar

# Save the color legend plot to PDF
dev.off()
