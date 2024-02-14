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
library("optparse")

args <- commandArgs(trailingOnly=T)

#args    
option_list = list(
  make_option(c("--GA"), type="character",
              help="GA", metavar="character"),
  make_option(c("--meta"), type="character",
              help="meta", metavar="character"),
  make_option(c("--OutputDir"), type="character",
              help="OutputDir", metavar="character"),
  make_option(c("--SampleName"), type="character",
              help="SampleName", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

WD <- opt$OutputDir
SampleName <- opt$SampleName
GA <- opt$GA
meta <- opt$meta

GBA <- read.delim(GA)
meta_data <- read.delim(meta)
GBA_filtered <- GBA[GBA$barcode %in% rownames(meta_data),]
GroupInfo <- data.frame(barcode=rownames(meta_data),
                        Ann=meta_data$Ann_v4)
head(GroupInfo)

GBA_filtered <- merge(x = GroupInfo, y = GBA_filtered,
                      by = "barcode", all.y = T)

dim(GBA_filtered)
head(gene_markers)
head(GBA_filtered)

GBA_filtered$Ann_gene_name <- paste0(GBA_filtered$Ann,"&",GBA_filtered$gene_name)
head(GBA_filtered)

SumTable <- GBA_filtered %>%
  group_by(Ann_gene_name) %>%
  summarise_at(vars(accessability), list(SumAcc = sum))
head(SumTable)

SumTable<- data.frame(celltype = as.character(lapply(strsplit(as.character(SumTable$Ann_gene_name),
                                                               split="&"), "[", 1)),
                       gene= as.character(lapply(strsplit(as.character(SumTable$Ann_gene_name),
                                                          split="&"), "[", 2)),
                       accessability = as.character(SumTable$SumAcc))
head(SumTable)
levels(factor(SumTable$celltype))
str(SumTable)
SumTable$accessability <- as.numeric(SumTable$accessability) 
SumTable_spread <- as_tibble(SumTable)
head(SumTable_spread)
setwd(WD)
SumTable_spread<- SumTable_spread %>% 
  spread(key = celltype, value = accessability,fill=0)
head(SumTable_spread)
dim(SumTable_spread)
write.table(SumTable_spread, file=paste0(SampleName,".GeneBodyACC.byGeneXCT.txt"), quote=F, row.names=F, col.names=T, sep="\t")

