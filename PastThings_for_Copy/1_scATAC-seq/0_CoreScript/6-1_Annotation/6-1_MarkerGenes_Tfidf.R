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

## 1) Calculate gene body accessibility
# load data
ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf"
chr <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"

bed_re1 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_Unique.bed"
#bed_re1 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/4_relk1_Unique.bed"

obj_re1 <- loadBEDandGenomeData(bed_re1, ann, chr)
GB_re1 <- GeneBodyAccessibility(obj_re1)

bed_re2 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_2_Unique.bed"
#bed_re2 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/4_relk1_2_Unique.bed"

obj_re2 <- loadBEDandGenomeData(bed_re2, ann, chr)
GB_re2 <- GeneBodyAccessibility(obj_re2)

head(GB_re2$sc_gene_ac)
dim(GB_re1$sc_gene_ac)

#Temp <- read.table(meta)
#tail(Temp)
#GBA_Re1 <- GB_re1$sc_gene_ac[which(GB_re1$sc_gene_ac$barcode %in% rownames(Temp)),]
#GBA_Re2 <- GB_re2$sc_gene_ac[which(GB_re2$sc_gene_ac$barcode %in% rownames(Temp)),]
#tail(GBA_Re2)
GBA_Re1 <- GB_re1$sc_gene_ac %>% filter(grepl('Zm', gene_name))
GBA_Re2 <- GB_re2$sc_gene_ac %>% filter(grepl('Zm', gene_name))
GBA <- rbind(GBA_Re1,GBA_Re2)
head(GBA)
tail(GBA)
dim(GBA)
write.table(GBA, file=paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_rel2.txt"), quote=F, row.names=F, col.names=T, sep="\t")
#write.table(GBA, file=paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_relk1.txt"), quote=F, row.names=F, col.names=T, sep="\t")

## MarkerGene Edit..
#Annotation <- read.table("/scratch/sb14489/3.scATAC/0.Data/MarkerGene/221110_EarMarker.txt", header=T)
#tail(Annotation)
#Bundle <- read.table("/scratch/sb14489/0.MarkerGene/Mesophyll_BundleSheath.txt", header=T)
#Bundle <- Bundle[which(Bundle$type=="bundle_sheath"),]
#New <- rbind(Annotation,Bundle)
#write.table(New, file=paste0("/scratch/sb14489/0.MarkerGene/Ear_InSitu.txt"), quote=F, row.names=F, col.names=T, sep="\t")


#### Start!
WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_whyzmCLE7DoesnotComeout"

Name <- "A619_ZmCLE7"

#Name <- "rel2AllMarkersWithCyc"
if (!dir.exists(WD)){
  dir.create(WD)
} else {
  print("Dir already exists!")
}
setwd(WD)

#### <- means recent ones
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"

geneact <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re.txt"

markers <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt"

pcs <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt" 
threads <- 10
target_cluster <- "LouvainClusters"
plot_each_CT<-"no"

dat <- loadData(meta, pcs, geneact, markers,target_cluster)
str(dat)
#dim(dat$marker.info)
head(dat$activity)[,c(1:10)]
#dat$activity["Zm00001eb000010",]
#dat$activity["Zm00001eb999999",]

#dat$marker.info <- dat$marker.info[c(1:9),]
#> head(dat$activity)[,c(1:10)]
#6 x 10 sparse Matrix of class "dgCMatrix"
#[[ suppressing 10 column names 'CB:Z:AAACGAAAGAAATCTG-2_rel2_2', 'CB:Z:AAACGAAAGGCAAGGG-2_rel2_2', 'CB:Z:AAACGAACAAGCGGTA-2_rel2_2' ... ]]
#Zm00001eb000010 . . . 34  . 12  .  .  . 12
#Zm00001eb000020 . . .  . 71  .  . 12 34 34
#Zm00001eb000050 . . .  .  .  .  .  . 12  .
#Zm00001eb000060 . . .  .  .  . 12  . 12  .
#Zm00001eb000070 . . .  .  .  .  .  . 12  .
#Zm00001eb000080 . . .  .  .  . 34  .  .  .
tail(dat$marker.info)
dim(dat$marker.info)
###############
##if we update the markers, we need to re-run this script to obtain the new set of smooth marker file
onlysmooth_marker <- "no" ##yes or no
run_denovo <- "no" ##yes or no
run_known <- "yes" ##yes or no

str(dat)
#saveRDS(dat,paste0(output_dir,"/GAobj.rds"))
#str(dat)
dim(dat$activity)
dim(dat$marker.info)
dim(dat$b)
dim(dat$h.pcs)

all.b <- dat$b
all.activity <- dat$activity
rownames(all.activity)[1:100]
all.hpcs <- dat$h.pcs
#dat$marker.info <- dat$marker.info[c(1:99),]
marker.info <- dat$marker.info
head(all.hpcs)
head(all.b)

message("The column name for meta data should be Cluster")
colnames(all.b)[length(colnames(all.b))] 
colnames(all.b)[length(colnames(all.b))] <- "Cluster"
dim(all.b)
head(all.b)
tail(marker.info)
dim(marker.info)

#tail(dat$marker.info )
#if (run_known == 'yes'){
dim(all.b)
dim(all.activity)
dim(all.hpcs)
dim(marker.info)

out <- runMajorPriori(all.b, all.activity , all.hpcs, marker.info, plot_each_CT,
                      marker_opt_dir,WD,onlysmooth_marker,threads=threads,output=Name)



