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


## Files
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,header=TRUE)

OutputDir <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4"

## Function1: Get nearest genes by cutoff of dACR.
Get_NearestGenes <-function(SelectedPeak){
  genes_ranges <- GRanges(seqnames=genes_data$V1, 
                        ranges=IRanges(start=genes_data$V2, end=genes_data$V3),
                        gene_id=genes_data$V4)

head(SelectedPeak)
split_names <- strsplit(SelectedPeak$Peak, split = "_")
chromosomes <- sapply(split_names, function(x) x[1])
starts <- as.integer(sapply(split_names, function(x) x[2]))
ends <- as.integer(sapply(split_names, function(x) x[3]))
granges_object <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))

nearest_genes <- distanceToNearest(granges_object, genes_ranges)
distances <- mcols(nearest_genes)$distance
nearest_genes_id <- genes_ranges[subjectHits(nearest_genes)]$gene_id
return(nearest_genes_id)
}
head(DEGInfo)
dim(DEGInfo)
NearestGenes <- Get_NearestGenes(DEGInfo)
length(NearestGenes)
## Save NearestGene Info!

GeneInfo <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.Curated.txt", fill = TRUE)
head(GeneInfo)
dim(GeneInfo)
ColNames <- GeneInfo[1,]
GeneInfo <- GeneInfo[-1, ]
colnames(GeneInfo) <- ColNames
head(GeneInfo)
head(GeneInfo$chr)
NearestGenes_df <- data.frame(gene_model = NearestGenes)
GeneInfo_unique <- GeneInfo[!duplicated(GeneInfo$gene_model), ]
NearestGenes_df$row_id <- seq_len(nrow(NearestGenes_df))
head(NearestGenes_df)
# Merge NearestGenes with GeneInfo
result <- merge(NearestGenes_df, GeneInfo_unique, by = "gene_model", all.x = TRUE)
# Reorder the merged result based on the row identifier
result <- result[order(result$row_id), ]
# Optionally, you can remove the row_id column if it's no longer needed
result$row_id <- NULL
# Check results
head(result)
DEGInfo_withGeneInfo <- cbind(DEGInfo,result)

write.table(DEGInfo_withGeneInfo, file=paste0(OutputDir,"/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion_NearestGENEINFO.txt"),
            quote=F, row.names=F, col.names=T, sep="\t")
