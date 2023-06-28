### 1) Cluster coloring again 
library(ggplot2)
library(RColorBrewer)
library(stringr)
head(loaded_meta_data)
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

library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)
library(GenomicRanges)

## R became wired..
BiocManager::install("cicero")
## Matrix update..
## rlang update
## recipes 
#* removing ‘/home/sb14489/.conda/envs/r_env/lib/R/library/recipes’
#* restoring previous ‘/home/sb14489/.conda/envs/r_env/lib/R/library/recipes’
#* removing ‘/home/sb14489/.conda/envs/r_env/lib/R/library/Seurat’
#* restoring previous ‘/home/sb14489/.conda/envs/r_env/lib/R/library/Seurat’

colorr <- colorRampPalette(brewer.pal(12,"Paired")[1:12])(12)
colorr <- colorr%>%str_replace_all(c("#4F96C4"="#060878",
                                     "#69BB54" ="#84f5d9",
                                     "#5D9E43"="#2d591b",
                                     "#DE9A89"="#e6d17e",
                                     "#EF595A"="#ce7ee6",
                                     "#FDA33F"="#876b58"))
colorr <- c(colorr,"#9d76d3","#734A12","#095119","#945151","#CA89AC")
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.REF_CELLs.metadata.txt"
loaded_meta_data <- read.table(meta)
head(loaded_meta_data)
levels(factor(loaded_meta_data$LouvainClusters))
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/MtHighCutoffAllMarker_RemoveBLSelected")

ggplot(loaded_meta_data, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("Ref_ColorRe_withSub.pdf", width=13, height=10)	


##################################
## 2) Subclustering ############
## Here I will use cluster 8
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref")
obj_A619_merged <- readRDS("Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.tfidf.rds")
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.REF_CELLs.metadata.txt"
loaded_meta_data <- read.table(meta)
head(loaded_meta_data)
dim(loaded_meta_data)
str(obj_A619_merged)
SelectedCluster <- loaded_meta_data[which(loaded_meta_data$LouvainClusters=="9"),]
dim(SelectedCluster)
str(obj_A619_merged)
obj_A619_merged$meta <- loaded_meta_data

## Remove cell cycle black list
blacklist_r <- read.table("/scratch/sb14489/3.scATAC/0.Data/CellCycle/CellCycle.txt")
head(blacklist_r)
blacklist_r <- blacklist_r[-1,]
blacklist.gr <- GRanges(seqnames=as.character(blacklist_r$V1),
                        ranges=IRanges(start=as.numeric(blacklist_r$V2),
                                       end=as.numeric(blacklist_r$V3)),
                        names=as.character(blacklist_r$V4))

head(blacklist.gr)
chr.seq.lengths_load <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai")
chr.seq.lengths <- as.numeric(chr.seq.lengths_load$V2)
names(chr.seq.lengths) <- chr.seq.lengths_load$V1
intervals <- tileGenome(seqlengths=chr.seq.lengths, tilewidth=500, cut.last.tile.in.chrom=TRUE)
head(intervals)
intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),] 
str(intervals)
regions <- as.data.frame(intervals)
head(regions)
regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
head(regions)

Temp <- obj_A619_merged$counts[rownames(obj_A619_merged$counts) %in% regions,]
dim(Temp)
dim(obj_A619_merged$counts)
obj_A619_merged$counts <- Temp

obj_A619_merged <- tfidf(obj_A619_merged, doL2=T)

Dir <- as.character("Ref") #NewDir name
Prefix <- "Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100"
WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/"
NumbeerOfWindow <- as.character(0)
SVDorNMF <-as.character("SVD")
NumberOfPC <- as.character(100)
NumbeerOfWindow <- as.character(0)

obj_A619_merged <- reduceDims(obj_A619_merged,method=SVDorNMF,
                              n.pcs=100,
                              cor.max=0.7,
                              num.var=as.numeric(NumbeerOfWindow),
                              verbose=T,
                              scaleVar=T,
                              doSTD=F,
                              doL1=F,
                              doL2=T,
                              refit_residuals=F)
#,svd_slotName="PCA"
#num.var=0



getwd()
#saveRDS(obj_A619_merged, file=paste0(out,".tfidf_ReducedD.rds"))
#obj_A619_merged <- readRDS(paste0(out,".tfidf_ReducedD.rds"))
###===== looks okay up to here ==========
str(obj_A619_merged)
# extract feature loadings and var/mean of tfidf
feat.data <- extractFeatureMetrics(obj_A619_merged, 0)
str(feat.data)
ids <- rownames(obj_A619_merged$PCA)

# remove batch effects with harmony --------------------------------------
dim(obj_A619_merged$PCA)
dim(obj_A619_merged$meta)
head(obj_A619_merged$meta)
head(obj_A619_merged$PCA)
obj_A619_merged$PCA <- obj_A619_merged$PCA[rownames(obj_A619_merged$meta),]
dim(obj_A619_mergedT)

ref.obj <- HarmonyMatrix(obj_A619_merged$PCA, meta_data=obj_A619_merged$meta, 
                         vars_use="library", do_pca=F,
                         #theta=c(3, 2), 
                         sigma=0.1, 
                         nclust=30,
                         max.iter.cluster=100,
                         max.iter.harmony=30,
                         return_object=T)


# create compressed harmony reference
getwd()
out <-  paste0("Ref_",Prefix,"_RemoveBLonlyMitoChloroChIP")

sym.ref <- symphony::buildReferenceFromHarmonyObj(ref.obj,
                                                  obj_A619_merged$meta,
                                                  feat.data$features,
                                                  feat.data$loadings,               
                                                  verbose = TRUE,         
                                                  do_umap = TRUE,       
                                                  umap_min_dist = 0.01,
                                                  save_uwot_path = paste0(out,'_uwot_model'))

obj_A619_merged$PCA <- t(ref.obj$Z_corr)
colnames(obj_A619_merged$PCA) <- paste0("PC_", 2:(ncol(obj_A619_merged$PCA)+1))
rownames(obj_A619_merged$PCA) <- rownames(obj_A619_merged$meta)
head(obj_A619_merged$PCA)
#saveRDS(sym.ref, file=paste0(out,".symphony.reference.rds"))
#saveRDS(obj_A619_merged, file=paste0(out,".afterHarmony.processed.rds"))
getwd()

## Here important!!

#sym.ref<- readRDS("Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.symphony.reference.rds")
#obj_A619_merged<- readRDS("Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.afterHarmony.processed.rds")


# reduce to 2-dimensions with UMAP ---------------------------------------
obj_A619_merged$UMAP <- sym.ref$umap$embedding
rownames(obj_A619_merged$UMAP) <- rownames(obj_A619_merged$PCA)
colnames(obj_A619_merged$UMAP) <- c("umap1","umap2")
str(obj_A619_merged)
#head(obj_A619_merged$counts)

str(obj_A619_merged)
SelectedCell <- rownames(SelectedCluster)
length(SelectedCell)
#head(obj_A619_merged$PCA_model)
#head(obj_A619_merged$residuals)
obj_sub <- list()
obj_sub$PCA <- obj_A619_merged$PCA[rownames(obj_A619_merged$PCA)%in%rownames(SelectedCluster),]
head(obj_sub$PCA)
obj_sub$UMAP <- obj_A619_merged$UMAP[rownames(obj_A619_merged$UMAP)%in%rownames(SelectedCluster),]
obj_sub$meta <- obj_A619_merged$meta[rownames(obj_A619_merged$meta)%in%rownames(SelectedCluster),]
obj_sub$counts <- obj_A619_merged$counts[,rownames(SelectedCluster)]
#head(obj_A619_merged$counts)[,c(1:10)]
dim(obj_sub$counts)
dim(obj_sub$meta)
# identify clusters using neighborhood graph -----------------------------
obj_A619_merged_Cluster <- callClusters(obj_sub, 
                                        res=0.2,
                                        k.near=3,
                                        verbose=T,
                                        min.reads=5e4,
                                        e.thresh=100,
                                        m.clst=25,
                                        cleanCluster=F,
                                        cl.method=2
                                        )
#m.clst=100
str(obj_A619_merged_Cluster)
pdf(paste0(out,"_Cluster9_Sub_res0.2_knear3.pdf"), width=10, height=10)
plotUMAP(obj_A619_merged_Cluster, cluster_slotName="Clusters", cex=0.2)
dev.off()
#dim(k)
str(obj_A619_merged_Cluster)
dim(obj_A619_merged_Cluster$Clusters)

obj_A619_merged_Cluster$Clusters$LouvainClusters <- paste0("Sub",obj_A619_merged_Cluster$Clusters$LouvainClusters)
head(obj_A619_merged_Cluster$Clusters)

Meta_Sub <- obj_A619_merged_Cluster$Clusters

head(loaded_meta_data)
rownames(Meta_Sub)
loaded_meta_data <- read.table(meta)
loaded_meta_data[rownames(Meta_Sub),]$LouvainClusters <- Meta_Sub$LouvainClusters
loaded_meta_data[c(20:40),]

levels(factor(loaded_meta_data$LouvainClusters))
loaded_meta_data <- loaded_meta_data[which(loaded_meta_data$LouvainClusters!="8"),]


############################################################
colorr <- colorRampPalette(brewer.pal(12,"Paired")[1:12])(12)
colorr <- colorr%>%str_replace_all(c("#4F96C4"="#060878",
                                     "#69BB54" ="#84f5d9",
                                     "#5D9E43"="#2d591b",
                                     "#DE9A89"="#e6d17e",
                                     "#EF595A"="#ce7ee6",
                                     "#FDA33F"="#876b58"))
colorr <- c(colorr,"#9d76d3","#734A12","#095119","#945151","#CA89AC")
head(loaded_meta_data)
loaded_meta_data <- loaded_meta_data[which(loaded_meta_data$LouvainClusters!="13"),]
loaded_meta_data <- loaded_meta_data[which(loaded_meta_data$LouvainClusters!="14"),]


levels(factor(loaded_meta_data$LouvainClusters))
Ann_V1 <- loaded_meta_data$LouvainClusters
Ann_V1 <- Ann_V1 %>%
  str_replace_all(c( 
    "11" = "SM_SPM_Base",
    "12" = "Meta_proto_phloem_sieve_element", 
    "SubSub1" = "Procambial_meristem", 
    "SubSub2" = "Procambial_meristem",
    "SubSub3" = "Meta_proto_xylem",
    "10"="Epidermis"))

Ann_V1 <- Ann_V1 %>%
  str_replace_all(c( 
    "1" = "Parenchyma",
    "2" = "Phloem_sieve_element",
    "3" = "Determinate_later_organs",
    "4" = "Sclerenchyma",
    "5" = "Bundle_Sheath",
    "6" = "Epidermis",
    "7" = "IM_SM_SPM",
    "9" = "Companion_cell"
    ))
Ann_V1 <- Ann_V1 %>%
  str_replace_all(c( "Epidermis" = "L1_Layer"))

loaded_meta_data$Ann_V1 <- Ann_V1
levels(factor(Ann_V1))

ggplot(loaded_meta_data, aes(x=umap1, y=umap2, color=factor(Ann_V1))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("Ref_ColorRe_withSub_withName.pdf", width=13, height=10)	
write.table(loaded_meta_data, file="Ref_withSub_withName.metadata.txt", quote=F, row.names=T, col.names=T, sep="\t")

