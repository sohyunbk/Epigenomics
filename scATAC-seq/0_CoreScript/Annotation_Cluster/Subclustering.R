library(ggplot2)
library(RColorBrewer)
library(stringr)
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
library("optparse")

args <- commandArgs(trailingOnly=T)
#args
option_list = list(
  make_option(c("--SampleName"), type="character",
              help="SampleName", metavar="character"),
  make_option(c("--MetaFile"), type="character",
              help="MetaFile", metavar="character"),
  make_option(c("--ObjAfterHarmony"), type="character",
              help="AfterHarmony.rds", metavar="character"),
  make_option(c("--AnnSlot"), type="character",
              help="AnnSlot", metavar="character"),
  make_option(c("--TargetClusterName"), type="character",
              help="TargetClusterName", metavar="character"),
  make_option(c("--OutputDir"), type="character",
              help="OutputDir", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
SampleName <- opt$SampleName
meta <- opt$MetaFile
obj_mergedFile <- opt$ObjAfterHarmony
Ann_Slot <- opt$AnnSlot
sCluster <- opt$TargetClusterName
OutputDir <- opt$OutputDir

## Ex - Bif3
#SampleName <- "bif3"
#setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/")
#meta <- "bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt"
#obj_merged <- readRDS("bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.rds")
#Ann_Slot <- "LouvainClusters"
#sCluster <- "1"
## Ex - rel2
#SampleName <- "rel2"
#setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/")
#meta <- "rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt"
#obj_merged <- readRDS("rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.rds")
#Ann_Slot <- "LouvainClusters"
#sCluster <- "1"

obj_merged <- readRDS(obj_mergedFile)
loaded_meta_data <- read.table(meta)
head(loaded_meta_data)
#loaded_meta_data[,Ann_Slot]

setwd(OutputDir)

### Function Start!! ###
RunSubClustering <- function(sCluster){
SelectedCluster <- loaded_meta_data[which(loaded_meta_data[,Ann_Slot]==sCluster),]
obj_sub <- list()
obj_sub$PCA <- obj_merged$PCA[rownames(SelectedCluster),]
obj_sub$UMAP <- obj_merged$UMAP[rownames(SelectedCluster),]
obj_sub$meta <- obj_merged$meta[rownames(SelectedCluster),]
obj_sub$counts <- obj_merged$counts[,rownames(SelectedCluster)]
#head(obj_A619_merged$counts)[,c(1:10)]
dim(obj_sub$counts)
dim(obj_sub$UMAP)
dim(obj_sub$PCA)
dim(obj_sub$meta)
colnames(obj_sub$UMAP) <- c("umap1","umap2")
# identify clusters using neighborhood graph -----------------------------
obj_Cluster <- callClusters(obj_sub,
                            res=0.2,
                            k.near=50,
                            verbose=T,
                            min.reads=5e4,
                            e.thresh=3,
                            m.clst=25,
                            cleanCluster=F,
                            cl.method=2)


ggplot(obj_Cluster$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.4) +
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(paste0("Cluster",sCluster,"_Sub_res0.2_knear50.pdf"), width=7, height=5)
write.table(obj_Cluster$Clusters, file=paste0("Cluster",sCluster,"_Sub_res0.2_knear50_Partmetadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")

## With different res and k.near option!
 obj_Cluster <- callClusters(obj_sub,
                            res=1,
                            k.near=100,
                            verbose=T,
                            min.reads=5e4,
                            e.thresh=3,
                            m.clst=25,
                            cleanCluster=F,
                            cl.method=2)
 ggplot(obj_Cluster$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
        geom_point(size=0.4) +
        theme_minimal()+
        guides(colour = guide_legend(override.aes = list(size=10)))
        ggsave(paste0("Cluster",sCluster,"_Sub_res1_knear100.pdf"), width=7, height=5)
 write.table(obj_Cluster$Clusters, file=paste0("Cluster",sCluster,"_Sub_res1_knear100_Partmetadata.txt"),
                        quote=F, row.names=T, col.names=T, sep="\t")
###### Make new UMAP <-- I am using the results from here!
obj_sub_reRun <- list()
obj_sub_reRun$meta <- obj_merged$meta[rownames(SelectedCluster),]
obj_sub_reRun$counts <- obj_merged$counts[,rownames(SelectedCluster)]
obj_sub_reRun$residuals <- obj_sub_reRun$counts
obj_sub_reRun$norm_method = "tfidf"

dim(obj_sub_reRun$counts)
SVDorNMF <-as.character("SVD")
obj_sub_reRun <- reduceDims(obj_sub_reRun,method=SVDorNMF,
                            n.pcs=100,
                            cor.max=0.7,
                            num.var=0,
                            verbose=T,
                            scaleVar=T,
                            doSTD=F,
                            doL1=F,
                            doL2=T,
                            refit_residuals=F)
obj_sub_reRun <- projectUMAP(obj_sub_reRun, verbose=T,
                             k.near=50, m.dist=0.01)
obj_sub_reRun <- callClusters(obj_sub_reRun,
                            res=1,
                            k.near=100,
                            verbose=T,
                            min.reads=5e4,
                            e.thresh=5,
                            m.clst=25,
                            cleanCluster=F,
                            cl.method=2)
str(obj_sub_reRun)
ggplot(obj_sub_reRun$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.4) +
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(paste0(SampleName,"_Cluster",sCluster,"_Recluster_Sub_res1_knear100.pdf"), width=7, height=5)


#pdf(paste0("Cluster",sCluster,"_Recluster_Sub_res0.2_knear50.pdf"), width=10, height=10)
#plotUMAP(obj_sub_reRun, cluster_slotName="Clusters", cex=0.2)
#dev.off()

head(obj_sub_reRun$meta)
str(obj_sub_reRun)
head(obj_sub_reRun$Clusters)
levels(factor(obj_sub_reRun$Clusters$LouvainClusters))
head(loaded_meta_data)
head(obj_sub_reRun$PCA)
Newmeta <- loaded_meta_data

obj_sub_reRun$Clusters$LouvainClusters <- paste0(sCluster,"_",obj_sub_reRun$Clusters$LouvainClusters)
Newmeta[rownames(obj_sub_reRun$Clusters),]$LouvainClusters <- obj_sub_reRun$Clusters$LouvainClusters
head(Newmeta)
write.table(Newmeta, file=paste0(SampleName,"_Cluster",sCluster,"_Recluster_Sub_res1_knear100_Allmetadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")
write.table(obj_sub_reRun$Clusters, file=paste0(SampleName,"_Cluster",sCluster,"_Recluster_Sub_res1_knear100_Partmetadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")
write.table(obj_sub_reRun$PCA, file=paste0(SampleName,"_Cluster",sCluster,"_Recluster_Sub_res1_knear100_Part.reduced_dimensions.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")

}

RunSubClustering(sCluster)


############################################################################
############################################################################
