library("here")
library(devtools)
library(Seurat)
library(stringr)
load_all('/home/sb14489/Socrates')
library(harmony)
library(symphony)
library(rlang)

setwd("/scratch/sb14489/3.scATAC_flo_Past/5.Socrates/4_CombineAll_AfterD")

## Load RDSs
ref_obj <- readRDS("Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.symphony.reference.rds")
obj_A619_merged <- readRDS("Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.afterHarmony.processed.rds")
m.obj <- readRDS("Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.ALL_CELLs_BeforeCluster.symphony_reference.rds")
kRes <- 0.7
kNear <- 30

str(ref_obj)
head(ref_obj$Z_corr)[,1]
head(m.obj$PCA)[,1]
str(m.obj)
head(m.obj$meta)

#row=sample col=PCs

# It seems like z_corr and PCs are same things.
tail(m.obj$PCA[which(m.obj$meta$SampleName=="A619"),])
Control_PCs <- m.obj$PCA[which(m.obj$meta$SampleName=="A619"),]
Query_PCs <- m.obj$PCA[which(m.obj$meta$SampleName!="A619"),]
dim(Control_PCs)
dim(Query_PCs)
table(m.obj$meta$SampleName)
#Control_PCs <- t(Control_PCs)
head(Control_PCs)
#Query_PCs <- t(Query_PCs)
obj_A619_merged$UMAP <- ref_obj$umap$embedding
rownames(obj_A619_merged$UMAP) <- rownames(obj_A619_merged$PCA)
colnames(obj_A619_merged$UMAP) <- c("umap1","umap2")

obj_A619_merged_Cluster <- callClusters(obj_A619_merged, 
                                        #res=0.3,
                                        res=0.7,
                                        #k.near=50,
                                        k.near=30,
                                        verbose=T,
                                        cleanCluster=F,
                                        cl.method=2,
                                        e.thresh=5)


## Mutants cluster by it self:
str(m.obj)
ex_obj <- m.obj
ex_obj2 <- 
ex_obj2$UMAP <- ex_obj$UMAP[which(ex_obj$meta$SampleName=="bif3"),]
ex_obj2$PCA <- ex_obj$PCA[which(ex_obj$meta$SampleName=="bif3"),]
ex_obj2$counts <- ex_obj$counts[,which(ex_obj$meta$SampleName=="bif3")]
colnames(ex_obj$counts)
rownames(ex_obj$PCA)
str(ex_obj)
ex_obj$meta <- ex_obj$meta[which(ex_obj$meta$SampleName=="bif3"),]

head(ex_obj$counts)
head(ex_obj$PCA)
str(obj_A619_merged)

Bif3_Cluster <- callClusters(ex_obj, 
                                  #res=0.3,
                                 res=0.7,
                                  #k.near=50,
                                 k.near=30,
                                verbose=T,
                              cleanCluster=F,
                            cl.method=2,
                             e.thresh=5)



tail(obj_A619_merged_Cluster$Clusters)
write.table(obj_A619_merged_Cluster$Clusters, file=paste0("Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.REF_CELLs.ResKN_",
kRes,"_",kNear,"_metadata.txt"), quote=F, row.names=F, col.names=T, sep="\t")

## Add here How to combine & adjust the clusters 

tail(Control_PCs)
levels(obj_A619_merged_Cluster$Clusters$LouvainClusters)
NewAnnotation <- obj_A619_merged_Cluster$Clusters$LouvainClusters
tail(NewAnnotation)
NewAnnotation <- sapply("Cluster", paste, NewAnnotation, sep="")
NewAnnotation <- sapply(NewAnnotation, paste, "Cluster", sep="")

NewAnnotation <- NewAnnotation %>%
str_replace_all(c( 
                    "Cluster11Cluster" = "Pith parenchyma",
                    "Cluster12Cluster" = "Xylem parenchyma", 
                  "Cluster13Cluster" = "Axillary meristem", 
                  "Cluster14Cluster" = "Phloem parenchyma", 
                  "Cluster15Cluster" = "Procambial meristem, Metaphloem, Protophloem", 
                  "Cluster16Cluster" = "Metaxylem, Protoxylem"))


NewAnnotation <- NewAnnotation%>%
  str_replace_all(c("Cluster1Cluster" = "NA1", 
                    "Cluster2Cluster" = "Phloem sieve element",
                    "Cluster3Cluster" = "Vascular Sclerenchyma", 
                    "Cluster4Cluster" = "Bundle Sheath", 
                    "Cluster5Cluster" = "SPM, SM",
                    "Cluster6Cluster" = "Floral meristem",
                    "Cluster7Cluster" = "L1 layer",
                    "Cluster8Cluster" = "NA2",
                    "Cluster9Cluster" = "IM, SPM, SM",
                    "Cluster10Cluster" = "L1 layer"))

levels(factor(NewAnnotation))


obj_A619_merged_Cluster$Clusters$LouvainClusters_Original <- obj_A619_merged_Cluster$Clusters$LouvainClusters
obj_A619_merged_Cluster$Clusters$LouvainClusters <- NewAnnotation
write.table(obj_A619_merged_Cluster$Clusters, 
            file=paste0("Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.REF_CELLs.ResKN_",
                kRes,"_",kNear,"_metadata_Annotation.txt"), quote=F, row.names=F, col.names=T, sep="\t")


## After clustering some samples in control are removed.
rownames(obj_A619_merged_Cluster$Clusters)
Control_PCs <- Control_PCs[rownames(obj_A619_merged_Cluster$Clusters),]


## Using the knn
# https://github.com/immunogenomics/symphony/blob/main/R/utils.R


knn_pred = class::knn(Control_PCs, Query_PCs, obj_A619_merged_Cluster$Clusters$LouvainClusters, k = 5, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
str(knn_pred)
str(m.obj)
head(m.obj$meta)
Clusters <- cbind(m.obj$meta,m.obj$UMAP)
head(Clusters)
head(obj_A619_merged_Cluster$Clusters)

Clusters_Control <- Clusters[rownames(obj_A619_merged_Cluster$Clusters),]
dim(Clusters_Control)
head(Clusters_Control)

Clusters_Control$Ref_Cluster <- obj_A619_merged_Cluster$Clusters$LouvainClusters
Clusters_Mutant <-Clusters[rownames(Query_PCs),]
Clusters_Mutant$Ref_Cluster <- knn_pred

head(Query_PCs)
head(Clusters_Mutant)
head(Clusters_Control)
MetaAll <- rbind(Clusters_Control,Clusters_Mutant)
head(MetaAll)

## Save Meta data
#write.table(MetaAll, file="Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.Clustering_BasedOnRefAllCells.metadata.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(MetaAll, file="Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.Clustering_BasedOnRefAllCells.Annotation.metadata.txt", quote=F, row.names=F, col.names=T, sep="\t")

## ggplot
library(ggplot2)
colorr <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(CluterTable$Ref_Cluster)))
colorr <- colorr%>%str_replace_all(c("#4F96C4"="#060878",
                                    "#69BB54" ="#84f5d9",
                                    "#5D9E43"="#2d591b",
                                    "#DE9A89"="#e6d17e",
                                    "#EF595A"="#ce7ee6",
                                    "#FDA33F"="#876b58"))



Clusters_Control$Ref_Cluster <- factor(Clusters_Control$Ref_Cluster,levels=c("Xylem parenchyma",
                                             "Pith parenchyma",
                                             "Metaxylem, Protoxylem",
                                             "Phloem parenchyma",
                                             "Procambial meristem, Metaphloem, Protophloem",
                                             "Phloem sieve element",
                                             "NA1",
                                             "NA2",
                                             "Bundle Sheath",
                                             "Vascular Sclerenchyma", 
                                             "Axillary meristem",
                                             "SPM, SM",
                                             "Floral meristem",
                                             "IM, SPM, SM",
                                             "L1 layer"))

ggplot(Clusters_Control, aes(x=umap1, y=umap2, color=Ref_Cluster)) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+
  guides(colour = guide_legend(override.aes = list(size=10)))+theme_minimal()

ggsave("Ref_clustering_withAnnotation.pdf", width=13, height=10)	

Clusters_Mutant$Ref_Cluster <- factor(Clusters_Mutant$Ref_Cluster,levels=c("Xylem parenchyma",
                                                                             "Pith parenchyma",
                                                                             "Metaxylem, Protoxylem",
                                                                             "Phloem parenchyma",
                                                                             "Procambial meristem, Metaphloem, Protophloem",
                                                                             "Phloem sieve element",
                                                                             "NA1",
                                                                             "NA2",
                                                                             "Bundle Sheath",
                                                                             "Vascular Sclerenchyma", 
                                                                             "Axillary meristem",
                                                                             "SPM, SM",
                                                                             "Floral meristem",
                                                                             "IM, SPM, SM",
                                                                             "L1 layer"))
ggplot(Clusters_Mutant, aes(x=umap1, y=umap2, color=Ref_Cluster)) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+
  guides(colour = guide_legend(override.aes = list(size=10)))+theme_minimal()

ggsave("Mutant_clustering_withAnnotation.pdf", width=13, height=10)	

MetaAll <- rbind(Clusters_Control,Clusters_Mutant)
ggplot(MetaAll, aes(x=umap1, y=umap2, color=Ref_Cluster)) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+
  guides(colour = guide_legend(override.aes = list(size=10)))+theme_minimal()

ggsave("AllCombined_clustering_withAnnotation.pdf", width=13, height=10)

CluterTable <- MetaAll
head(CluterTable)
levels(CluterTable$Ref_Cluster)
#x <- c(1,2,3,6,3,6)
#y <- c(4,5,7,8,9,0)
#dist(cbind(x,y))

i <- "6"
table(CluterTable$SampleName)[["A619"]]
head(CluterTable)
#colorr <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(CluterTable$Ref_Cluster)))

Plotlist <- list()
ggplotData <- list()
i <- levels(CluterTable$Ref_Cluster)[1]
k = 1
for (i in levels(CluterTable$Ref_Cluster)){
  
  SubCluster <- CluterTable[which(CluterTable$Ref_Cluster==i),]
  #head(SubCluster)
  Plotlist[[i]] <- data.frame()
  for (j in levels(factor(SubCluster$SampleName))){
    Count <- nrow(SubCluster[which(SubCluster$SampleName==j),])
    Temp <- data.frame(SampleName=j,Fre=Count,
                       Ratio=(Count/table(CluterTable$SampleName)[[j]])*100)
    #rownames(Temp) <- NULL
    Plotlist[[i]] <- rbind(Plotlist[[i]],Temp)
  }
  Plotlist[[i]]$SampleName <- factor(Plotlist[[i]]$SampleName ,levels=c("relk1","rel2","bif3","A619"))
  
  ggplotData[[i]] <- ggplot(data=Plotlist[[i]], aes(x=SampleName, y=Ratio)) +theme_minimal()+
    geom_bar(stat="identity",fill=colorr[k])+coord_flip()+
    xlab("") +
    ylab("Ratio(%)")+ggtitle(i)+ theme(plot.title = element_text(size = 10))
  k=k+1
}
library(easyGgplot2)

tiff("Ratio_by_Cluster_Barplot.tiff", width = 10, height = 7, units = 'in', res = 300)
ggplot2.multiplot(ggplotData[[1]], ggplotData[[2]], ggplotData[[3]], ggplotData[[4]],
       ggplotData[[5]], ggplotData[[6]], ggplotData[[7]], ggplotData[[8]],
          ggplotData[[9]], ggplotData[[10]], ggplotData[[11]], ggplotData[[12]],
       ggplotData[[13]], ggplotData[[14]], ggplotData[[15]],
                cols=4)
dev.off()


######################################
CluterTable <- m.obj$Clusters

head(DistanceMatrix)

#x <- c(1,2,3,6,3,6)
#y <- c(4,5,7,8,9,0)
#dist(cbind(x,y))

i <- "6"
table(CluterTable$SampleName)[["A619"]]

Plotlist <- list()
ggplotData <- list()

library(ggplot2)
Colorr <- c ("#E41A1C" ,"#FF7F00", "#7FC97F"  ,"#377EB8", "#984EA3","#F781BF")

for (i in levels(CluterTable$LouvainClusters)){
  SubCluster <- CluterTable[which(CluterTable$LouvainClusters==i),]
  #head(SubCluster)
  Plotlist[[i]] <- data.frame()
  for (j in levels(factor(SubCluster$SampleName))){
    Count <- nrow(SubCluster[which(SubCluster$SampleName==j),])
    Temp <- data.frame(SampleName=j,Fre=Count,
                       Ratio=(Count/table(CluterTable$SampleName)[[j]])*100)
    #rownames(Temp) <- NULL
    Plotlist[[i]] <- rbind(Plotlist[[i]],Temp)
  }
  ggplotData[[i]] <- ggplot(data=Plotlist[[i]], aes(x=SampleName, y=Ratio)) +
    geom_bar(stat="identity",fill=)+coord_flip()
}




## Outputtable:
# Cluster# <- list Sample# Number# Ratio 
