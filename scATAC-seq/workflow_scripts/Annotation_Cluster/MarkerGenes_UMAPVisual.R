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
source("/home/sb14489/Epigenomics/scATAC-seq/0_Function/GeneBodyAccessibility.R")
source("/home/sb14489/Epigenomics/scATAC-seq/0_Function/MarkerGenes_Tfidf_functions.R")
library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)
##################################################################################
#### Start!
library("optparse")
library(rlang)
library(ggplot2)
option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name"),
  make_option(c("--meta"), type="character",
              help="meta", metavar="character"),
  make_option(c("--geneact"), type="character",
              help="geneact", metavar="character"),
  make_option(c("--markers"), type="character",
              help="markers", metavar="character"),
  make_option(c("--pcs"), type="character",
              help="pcs", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

WD <- opt$WD
Name <- opt$Name
meta <- opt$meta
geneact <- opt$geneact
markers <- opt$markers
pcs <- opt$pcs

if (!dir.exists(WD)){
  dir.create(WD)
} else {
  print("Dir already exists!")
}
setwd(WD)

threads <- 20
target_cluster <- "LouvainClusters"
plot_each_CT<-"no"

dat <- loadData(meta, pcs, geneact, markers, target_cluster)
str(dat)
dat$marker.info <- dat$marker.info
###############
##if we update the markers, we need to re-run this script to obtain the new set of smooth marker file
onlysmooth_marker <- "no" ##yes or no
run_denovo <- "no" ##yes or no

str(dat)

all.b <- dat$b
all.activity <- dat$activity
#rownames(all.activity)[1:100]
all.hpcs <- dat$h.pcs
marker.info <- dat$marker.info

head(all.b)
message("The column name for meta data should be Cluster")
colnames(all.b)[length(colnames(all.b))]
colnames(all.b)[length(colnames(all.b))] <- "Cluster"
head(all.b)
tail(dat$marker.info)

out <- runMajorPriori(all.b, all.activity , all.hpcs, marker.info, plot_each_CT,
                      marker_opt_dir,WD,onlysmooth_marker,threads=threads,output=Name)
saveRDS(out,paste0(WD,'/MarkerGenes.rds'))
