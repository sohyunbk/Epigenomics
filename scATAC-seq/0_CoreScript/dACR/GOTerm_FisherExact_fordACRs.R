library(devtools)
library("fgsea")
library("here")
library(devtools)
library(tidyverse)
library(Matrix)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library("optparse")
library(GenomicRanges)
library(topGO)
library(ggplot2)


option_list = list(
  make_option(c("--GOData"), type="character",
              help="GOData", metavar="character"),
  make_option(c("--WDir_fordACR"), type="character",
              help="WDir_fordACR", metavar="character"),
  make_option(c("--FileNameFix"), type="character",
              help="FileNameFix", metavar="character"),
  make_option(c("--OutPutDir"), type="character",
              help="OutPutDir", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Start here 
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
#acr_data <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.FDR0.01Bif3Higher.bed")
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,header=TRUE)
GOData <- "/scratch/sb14489/0.Reference/Maize_B73/GOTerm/B73_GO.RemoveBlank.out"

MakeGOlist <-function(GOTerm){
  zm_func_genes <- read.table(GOData,
                              header=TRUE)
  BP <- zm_func_genes[which(zm_func_genes$ontology==GOTerm),]
  #head(BP)
  #dim(BP)
  BP$gene_id <- sapply(strsplit(as.character(BP$qpid),'_'), "[", 1) 
  BP$Name <- paste0(BP$desc,".",BP$goid)
  BP_list <- list()
  BPNames <- levels(factor(BP$Name))
  length(BPNames)
  i <- BPNames[2]
  for(i in BPNames){
    Sub <- BP[which(BP$Name == i),]
    #head(Sub)
    BP_list[[i]] <- Sub$gene_id
  }
  return(BP_list)}


# Convert to GRanges

genes_ranges <- GRanges(seqnames=genes_data$V1, 
                        ranges=IRanges(start=genes_data$V2, end=genes_data$V3),
                        gene_id=genes_data$V4)

#GOList_BP<- MakeGOlist("BP")
head(DEGInfo)
SelectedPeak <- DEGInfo[DEGInfo$FDR < 0.05 & DEGInfo$logFC > 0,]
head(SelectedPeak)
split_names <- strsplit(SelectedPeak$Peak, split = "_")
chromosomes <- sapply(split_names, function(x) x[1])
starts <- as.integer(sapply(split_names, function(x) x[2]))
ends <- as.integer(sapply(split_names, function(x) x[3]))
granges_object <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))

nearest_genes <- distanceToNearest(granges_object, genes_ranges)
distances <- mcols(nearest_genes)$distance
nearest_genes_id <- genes_ranges[subjectHits(nearest_genes)]$gene_id
head(nearest_genes_id)

FisherExactTest <-function(GOList_BP){
Pvalues <-c()
GOName <- c()
for (i in c(1:length(GOList_BP))){
  
##########
## Fisher Exact Test

#   DEG   NonDEG
#GO 2   63
#NonGO 1437  
CheckedGO <- GOList_BP[[i]]  
GONames <- names(GOList_BP)[i]
## it has repetitive number! 
CheckedGO_unique <- unique(CheckedGO)

GO_DEG <- length(intersect(CheckedGO_unique, nearest_genes_id))
NonGO_DEG <- length(nearest_genes_id)-GO_DEG

NonDEGs <- setdiff(genes_ranges$gene_id, nearest_genes_id)

GO_NonDEG <- length(intersect(CheckedGO_unique, NonDEGs))
NonGO_NonDEG <- length(genes_ranges$gene_id)-GO_NonDEG

dat <- data.frame(
  "DEG" = c(GO_DEG, NonGO_DEG),
  "NonDEG" = c(GO_NonDEG, NonGO_NonDEG),
  row.names = c("GO", "NonGO"),
  stringsAsFactors = FALSE
)

Pvalues <- c(Pvalues,fisher.test(dat)$p.value)
GOName <-c(GOName,GONames)
}

FDR<-p.adjust(Pvalues, method="BH")
Result <- data.frame(GOName=GOName,Pvalues=Pvalues,FDR=FDR)
Sig <- Result[Result$FDR <0.05,]
Sig <- Result[Result$Pvalues <0.001,]
Sig_ordered <- Sig[rev(order(Sig$Pvalues)), ]
Sig_ordered$logP <- -log10(Sig_ordered$Pvalues)

ggplot(Sig_ordered, aes(y = reorder(GOName, logP), x = logP)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  labs(y = "GO Term", x = "-log10(P-value)", title = "BP_IM-OC_Bif3Higher")

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/")
ggsave("GO_Bif3Higher_Pvalue0.001.pdf", width = 10, height = 2, dpi = 300)
