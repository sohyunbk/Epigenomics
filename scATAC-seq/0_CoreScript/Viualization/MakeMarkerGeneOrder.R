# Assuming you have already installed the necessary packages
library(ggplot2)

# Sample data (replace with your actual data)
set.seed(123) # For reproducibility
## In RNA-seq, the dot size is persent expression and Color is average expression.
## In ATAC-seq, the dot size is persent acc and Color is  average acc or zscore.
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
library(stringr)
library("optparse")

option_list = list(
  make_option(c("--MarkerGene"), type="character",
              help="MarkerGene", metavar="character"),
  make_option(c("--OutPathandPrefix"), type="character",
              help="OutPathandPrefix", metavar="character"),
  make_option(c("--CellOrdertxt"), type="character",
              help="CellOrdertxt", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
gene <- opt$MarkerGene
OutputPathandName<- opt$OutPathandPrefix
CellOrder <- opt$CellOrdertxt

#gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt"
#OutputPathandName <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_GeneNameOrder"
#CellOrder <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt"


gene_markers <- read.delim(gene)
gene_markers <- gene_markers  %>%
  arrange(type)
##########
head(gene_markers)
CellTypeOrder <- readLines(CellOrder)
# Perform a left join to add the "name" column to the combined_table based on geneID

first_match <- function(name, order_list) {
  for (order_item in order_list) {
    if (str_detect(name, order_item)) {
      return(order_item)
    }
  }
  return(NA)  # Return NA if no match is found
}
# Create a new column for sorting 
## type is cell type matched with the marker gene!!
## Ann is Annotated cell type
gene_markers$sorting_key <- sapply(gene_markers$type, first_match, order_list = CellTypeOrder)
head(gene_markers)
# Convert the sorting_key column to a factor with levels in the order of CellTypeOrder
gene_markers$sorting_key <- factor(gene_markers$sorting_key, levels = CellTypeOrder)
# Arrange the table by the sorting_key column
gene_markers <- gene_markers %>% 
  arrange(sorting_key, type)
head(gene_markers)

OrderedName <- unique(gene_markers$name)
writeLines(OrderedName, paste0(OutputPathandName,".txt"))
