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

GettingGOResult <-function(BP_list,SelectedPeak,GOTerm){
  #str(BP_list)
  #CellName <- "XylemParenchyma"
  
  selected <- SelectedPeak  %>% 
    dplyr::select(Peak, FDR)
  #head(zm_bs_denovo.selected)
  ranks <- deframe(selected)
  ranks_ordered <- rev(ranks[order(ranks)])
  Peak <- names(ranks_ordered)
  # Split the names into components
  split_names <- strsplit(Peak, split = "_")
  
  # Extract chromosomes, start and end positions
  chromosomes <- sapply(split_names, function(x) x[1])
  starts <- as.integer(sapply(split_names, function(x) x[2]))
  ends <- as.integer(sapply(split_names, function(x) x[3]))
  granges_object <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))
  
  nearest_genes <- distanceToNearest(granges_object, genes_ranges)
  distances <- mcols(nearest_genes)$distance
  nearest_genes_id <- genes_ranges[subjectHits(nearest_genes)]$gene_id
  
  names(ranks_ordered) <- nearest_genes_id
  
  #head(ranks_ordered)
  #length(ranks_ordered)
  # probalby sort ranks
  
  # stats : Named vector of gene-level stats.
 
  #https://code.bioconductor.org/browse/fgsea/RELEASE_3_18/
  fgseaRes <- fgsea(pathways=BP_list, stats=ranks_ordered, 
                    nperm=10000,
                    minSize = 10, maxSize = 600)
  #str(fgseaRes)  
  #head(fgseaRes)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(NES)
  
  Output <- fgseaResTidy[order(fgseaResTidy$pval,decreasing=FALSE,na.last=TRUE),]
  #head(Output)
  OutputTable <- data.frame(Output)
  #colnames(OutputTable)
  OutputTable <- subset(OutputTable, select = -c(leadingEdge))
  #str(OutputTable)
  #head(OutputTable)
  setwd(opt$OutPutDir)
  CellName <-  sub(FileEnd, "", FileName)
  write.table(OutputTable, file=paste0(CellName,".",GOTerm,".txt"),
              quote=F, row.names=F, col.names=T, sep="\t")
}


## Start here 
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
acr_data <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.FDR0.01Bif3Higher.bed")
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")

# Convert to GRanges
acr_ranges <- GRanges(seqnames=acr_data$V1, 
                      ranges=IRanges(start=acr_data$V2, end=acr_data$V3))

genes_ranges <- GRanges(seqnames=genes_data$V1, 
                        ranges=IRanges(start=genes_data$V2, end=genes_data$V3),
                        gene_id=genes_data$V4)

# Find distance to nearest gene
nearest_genes <- distanceToNearest(acr_ranges, genes_ranges)
distances <- mcols(nearest_genes)$distance
nearest_genes_id <- genes_ranges[subjectHits(nearest_genes)]$gene_id

GOData <- "/scratch/sb14489/0.Reference/Maize_B73/GOTerm/B73_GO.RemoveBlank.out"
GOList_BP<- MakeGOlist("BP")

DEGInfo <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion.txt",header=TRUE)
head(DEGInfo)
SelectedPeak <- DEGInfo[DEGInfo$FDR < 0.05 & DEGInfo$logFC >0 ,]
head(SelectedPeak)

SelectedPeak$Peak

for (File in files){
  GOList_BP<- MakeGOlist("BP")
  GettingGOResult(GOList_BP,File,"BP")
  GOList_CC<- MakeGOlist("CC")
  GettingGOResult(GOList_CC,File,"CC")
  GOList_MF<- MakeGOlist("MF")
  GettingGOResult(GOList_MF,File,"MF")
}
