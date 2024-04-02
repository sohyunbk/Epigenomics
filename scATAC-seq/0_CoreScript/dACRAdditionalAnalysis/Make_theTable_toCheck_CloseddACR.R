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
library(edgeR)
library(preprocessCore)
library(GenomicRanges)
library(gplots)

#### 1) Find the two closest genes!
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion_NearestGENEINFO.txt"
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,fill=TRUE,header=TRUE)
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC")
FimoWD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05A619Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/"

Fimo1 <- paste0(FimoWD,"/fimo_out_1/fimo.tsv")
Fimo3 <- paste0(FimoWD,"/fimo_out_3/fimo.tsv")
Fimo4 <- paste0(FimoWD,"/fimo_out_4/fimo.tsv")
Fimo5 <- paste0(FimoWD,"/fimo_out_5/fimo.tsv")
Fimo7 <- paste0(FimoWD,"/fimo_out_7/fimo.tsv")


Make_heatmap(Fimo1,"Fimo1_GCACAGCAGCR")
Make_heatmap(Fimo3,"Fimo3_GCAGCATGC")
Make_heatmap(Fimo4,"Fimo4_CGCGCCGCGCC")
Make_heatmap(Fimo5,"Fimo5_YAGAGAGAGA")
Make_heatmap(Fimo7,"Fimo7_GCTAGCTAGC")
FimoFile <- Fimo1
Make_heatmap <- function(FimoFile,OutputName){

DEGInfo_A619Higher <- DEGInfo[(DEGInfo$FDR < 0.05) & (DEGInfo$logFC < 0),]
FimoOut <- read.table(FimoFile,header=TRUE)
## Find overlap.
DEGInfo$chr <- sapply(strsplit(as.character(DEGInfo$Peak), "_"), `[`, 1)
DEGInfo$start <- as.integer(sapply(strsplit(as.character(DEGInfo$Peak), "_"), `[`, 2))
DEGInfo$end <- as.integer(sapply(strsplit(as.character(DEGInfo$Peak), "_"), `[`, 3))
DEG_ranges <- GRanges(
  seqnames = DEGInfo$chr,
  ranges = IRanges(start = DEGInfo$start, end = DEGInfo$end)
)

FimoOut_ranges <- GRanges(
  seqnames = FimoOut$sequence_name,
  ranges = IRanges(start = FimoOut$start, end = FimoOut$stop)
)
overlaps <- findOverlaps(DEG_ranges, FimoOut_ranges)
overlapping_DEGInfo <- DEGInfo_A619Higher[queryHits(overlaps), ]
#overlapping_FimoOut <- FimoOut[subjectHits(overlaps), ]
dim(overlapping_DEGInfo)
head(overlapping_DEGInfo)
overlapping_DEGInfo <- unique(overlapping_DEGInfo)

## load gene*cell table.
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
##############
SelectedGenes <- overlapping_DEGInfo$gene_model

CTOrder <- gsub("-", ".", CTOrder)
new_CTOrder <- unlist(lapply(CTOrder, function(x) c(paste("A619&", x, sep=""), paste("Bif3&", x, sep=""))))
dim(qnorm_data)
qnorm_data_orderd <- qnorm_data[,new_CTOrder]
QNorm_SelectedGene <- qnorm_data_orderd[rownames(qnorm_data_orderd)%in%SelectedGenes,]
dim(QNorm_SelectedGene)
head(QNorm_SelectedGene)

# Define color palette for the heatmap
my_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(255)
# Create heatmap


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
### box plot ##########
FCTable_df <- as.data.frame(FCTable_ordered)
CTOrder <- gsub("\\.", "_", CTOrder)
colnames(FCTable_df) <- factor(colnames(FCTable_df),levels=CTOrder)
head(FCTable_df)
cols_to_keep <- !grepl("Unknown|G2", colnames(FCTable_df))
FCTable_df_NotUnknown <- FCTable_df[, cols_to_keep]
long_data <- gather(FCTable_df_NotUnknown, key = "CellType", value = "Expression")
long_data$CellType <- factor(long_data$CellType,levels=CTOrder[cols_to_keep])
# Create the box plot
dim(long_data)
head(long_data)
getwd()
colorr <- c("#4F96C4","#84f5d9","#0bd43d","#d62744","#FDA33F","#060878","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#adafde","#DE9A89","#5703ff",
            "#deadce","#fc53b6")
ggplot(long_data, aes(x = CellType, y = Expression, fill = CellType)) +
  geom_boxplot(size = 0.1) +  # Adjust size here for thinner lines
  scale_fill_manual(values = colorr) +
  labs(y = "log FC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(color = NULL))) +
  scale_y_continuous(limits = c(-2, 2))
  
ggsave(paste0("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_",OutputName,".pdf"), width=15, height=15)

write.table(FCTable_ordered,paste0("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_",OutputName,".txt"), 
            quote=F, row.names=T, col.names=T, sep="\t")

}

