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
library(edgeR)
library(parallel)
library(png)
library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)

Get_Peak_Tn5Sparse_NonBinary <- function(bed,Peak){
    Tn5_bed <- read.table(bed)
    Peak_bed <- read.table(Peak)
    ## Load Tn5 file with GRanges
    Tn5_Grange <-  GRanges(seqnames = Tn5_bed$V1,
                           ranges = IRanges(start = Tn5_bed$V2,
                           end = Tn5_bed$V3,
                           names = Tn5_bed$V4))
    Peak_bed$Name <- paste(Peak_bed$V1,Peak_bed$V2,Peak_bed$V3,sep="_")
    Peak_Grange <-  GRanges(seqnames = Peak_bed$V1,
                           ranges = IRanges(start = Peak_bed$V2,
                           end = Peak_bed$V3,
                           names = Peak_bed$Name))    
    
    ## Find overlap
    hits_Within <- findOverlaps(Tn5_Grange,  Peak_Grange,minoverlap=1,
                                type=c("within"),select="all",ignore.strand = TRUE)
    
    Intersect <- paste(names(Peak_Grange)[hits_Within@to],
                       names(Tn5_Grange)[hits_Within@from],sep="/_Com_/")
    
    Intersect <- table(Intersect)
    Intersect<- data.frame(gene_name = as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                 split="/_Com_/"), "[", 1)),
                        barcode= as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                              split="/_Com_/"), "[", 2)),
                        accessability = as.character(Intersect))
    
      return(Intersect)}

##################################################################################


library("optparse")
library(rlang)
library(ggplot2)
option_list = list(
  make_option(c("--Peak"), type="character",
              help="Peak bed format", metavar="character"),
  make_option(c("--Re1_bed"), type="character",
              help="Re1_bed", metavar="character"),
  make_option(c("--Re2_bed"), type="character",
              help="Re2_bed", metavar="character"),
  make_option(c("--OutFileName"), type="character",
              help="OutFileName", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
## We need two beds as we have two replicates! :)
Peak <- opt$Peak
#chr <- opt$ChrFai
bed_re1 <- opt$Re1_bed
bed_re2 <- opt$Re2_bed
OutfileName <- opt$OutFileName

# load data
#Peak <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed"
#chr <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"
#bed_re1 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_Unique.bed"
#bed_re2 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_2_Unique.bed"
#bed <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/Test.bed"

Peak_Barcode_Tn5_Re1 <- Get_Peak_Tn5Sparse_NonBinary(bed_re1,Peak)
Peak_Barcode_Tn5_Re2 <- Get_Peak_Tn5Sparse_NonBinary(bed_re2,Peak)

Peak_Barcode_Tn5 <- rbind(Peak_Barcode_Tn5_Re1,Peak_Barcode_Tn5_Re2)
write.table(Peak_Barcode_Tn5, file=OutfileName, quote=F, row.names=F, col.names=T, sep="\t")

