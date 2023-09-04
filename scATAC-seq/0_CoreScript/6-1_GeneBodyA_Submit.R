library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
#library(varistran)
library(edgeR)
library(parallel)
library(png)
source("/home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/0_Functions/GeneBodyAccessibility.R")
source("/home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/0_Functions/MarkerGenes_Tfidf_functions.R")
library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)

##################################################################################
## Keep original way to do it but naming should be changed
#conda activate r_env
#WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/MarkerGene"

## 1) Calculates gene body accessibility
library("optparse")
library(rlang)
library(ggplot2)
option_list = list(
  make_option(c("--Ann"), type="character",
              help="Ann_gtf", metavar="character"),
  make_option(c("--ChrFai"), type="character",
              help="ChrFai"),
  make_option(c("--Re1_bed"), type="character",
              help="Re1_bed", metavar="character"),
  make_option(c("--Re2_bed"), type="character",
              help="Re2_bed", metavar="character"),
  make_option(c("--OutFileName"), type="character",
              help="OutFileName", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ann <- opt$Ann
chr <- opt$ChrFai
bed_re1 <- opt$Re1_bed
bed_re2 <- opt$Re2_bed
OutfileName <- opt$OutFileName

# load data
#ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf_AddZmCLE7.gtf"
#chr <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"
#bed_re1 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_Unique.bed"
#bed_re1 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/4_relk1_Unique.bed"
#bed_re2 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_2_Unique.bed"
#bed_re2 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/4_relk1_2_Unique.bed"
#
obj_re1 <- loadBEDandGenomeData(bed_re1, ann, chr)
GB_re1 <- GeneBodyAccessibility(obj_re1)

obj_re2 <- loadBEDandGenomeData(bed_re2, ann, chr)
GB_re2 <- GeneBodyAccessibility(obj_re2)
#subset(transcriptLengths(obj_re1$gff, with.cds_len = TRUE), gene_id == "Zm00001eb999999")
head(GB_re2$sc_gene_ac)
dim(GB_re1$sc_gene_ac)

GBA_Re1 <- GB_re1$sc_gene_ac %>% filter(grepl('Zm', gene_name))
GBA_Re2 <- GB_re2$sc_gene_ac %>% filter(grepl('Zm', gene_name))
GBA <- rbind(GBA_Re1,GBA_Re2)
head(GBA)
tail(GBA)
dim(GBA)
write.table(GBA, file=OutfileName, quote=F, row.names=F, col.names=T, sep="\t")
