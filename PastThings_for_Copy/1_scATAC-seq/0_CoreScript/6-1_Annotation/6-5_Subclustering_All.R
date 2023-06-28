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


### 1) Cluster coloring again 
colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#de62b9","#0bd43d","#cf1d35","#deadce","#adafde")

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/")
meta <- "Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.metadata.txt"
loaded_meta_data <- read.table(meta)
head(loaded_meta_data)
levels(factor(loaded_meta_data$LouvainClusters))

ggplot(loaded_meta_data, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("Ref_Color_ggplot.pdf", width=13, height=10)	


##################################
## 2) Subclustering #0###########
## Here I will use cluster 8

obj_A619_merged <- readRDS("Ref_RemoveBLonlyMitoChloroChIP.ALL_CELLs_BeforeCluster.symphony_reference.rds")
sym.ref <- readRDS("Ref_RemoveBLonlyMitoChloroChIP.symphony.reference.rds")

head(loaded_meta_data)
str(obj_A619_merged)
head(obj_A619_merged$counts)[,c(1:10)]
head(obj_A619_merged$UMAP)
head(sym.ref$umap$embedding)
obj_A619_merged$UMAP <- sym.ref$umap$embedding

## DefineFunction: 
sCluster <- "1"
RunSubClustering <- function(sCluster){
SelectedCluster <- loaded_meta_data[which(loaded_meta_data$LouvainClusters==sCluster),]
obj_sub <- list()
obj_sub$PCA <- obj_A619_merged$PCA[rownames(SelectedCluster),]
obj_sub$UMAP <- obj_A619_merged$UMAP[rownames(SelectedCluster),]
obj_sub$meta <- obj_A619_merged$meta[rownames(SelectedCluster),]
obj_sub$counts <- obj_A619_merged$counts[,rownames(SelectedCluster)]
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

###### Make new UMAP
obj_sub_reRun <- list()
obj_sub_reRun$meta <- obj_A619_merged$meta[rownames(SelectedCluster),]
obj_sub_reRun$counts <- obj_A619_merged$counts[,rownames(SelectedCluster)]
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
ggsave(paste0("Cluster",sCluster,"_Recluster_Sub_res1_knear100.pdf"), width=7, height=5)	


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
write.table(Newmeta, file=paste0("Cluster",sCluster,"_Recluster_Sub_res1_knear100_Allmetadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")
write.table(obj_sub_reRun$Clusters, file=paste0("Cluster",sCluster,"_Recluster_Sub_res1_knear100_Partmetadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")
write.table(obj_sub_reRun$PCA, file=paste0("Cluster",sCluster,"_Recluster_Sub_res1_knear100_Part.reduced_dimensions.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")

}
RunSubClustering("1")
RunSubClustering("3")
RunSubClustering("4")
RunSubClustering("7")

############################################################################
############################################################################







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

