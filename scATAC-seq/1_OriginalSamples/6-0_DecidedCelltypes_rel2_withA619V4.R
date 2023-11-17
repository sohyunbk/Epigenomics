## Make metafile & plot for estimated annotation
library(ggplot2)
library(stringr)

## Make new metafile bycombining all the meta
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/rel2/")
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt"
loaded_meta_data <- read.table(meta)
head(loaded_meta_data)

#For rel2, I did sub clustering for cluster 1,3,4
SubClusterDir <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/rel2/"
Cluster1 <- read.table(paste0(SubClusterDir,"Cluster1_Sub_res1_knear100_Partmetadata.txt"))
Cluster1$LouvainClusters <- paste0("1_",Cluster1$LouvainClusters)
Cluster3 <- read.table(paste0(SubClusterDir,"Cluster3_Sub_res1_knear100_Partmetadata.txt"))
Cluster3$LouvainClusters <- paste0("3_",Cluster3$LouvainClusters)
Cluster4 <- read.table(paste0(SubClusterDir,"Cluster4_Sub_res1_knear100_Partmetadata.txt"))
Cluster4$LouvainClusters <- paste0("4_",Cluster4$LouvainClusters)
head(Cluster1)
#levels(factor(loaded_meta_data$LouvainClusters))
loaded_meta_data[rownames(Cluster1),]$LouvainClusters <- Cluster1$LouvainClusters
loaded_meta_data[rownames(Cluster3),]$LouvainClusters <- Cluster3$LouvainClusters
loaded_meta_data[rownames(Cluster4),]$LouvainClusters <- Cluster4$LouvainClusters

levels(factor(loaded_meta_data$LouvainClusters))

NewMeta <- loaded_meta_data
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="1"),]
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="3"),]
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="4"),]

dim(loaded_meta_data)
dim(NewMeta)
levels(factor(NewMeta$LouvainClusters))

NewMeta$Ann_v4 <- "Ann_v4"
levels(factor(NewMeta$Ann_v4))

NewMeta[which(NewMeta$LouvainClusters =="1_1"),]$Ann_v4 <- "Unknown1"
NewMeta[which(NewMeta$LouvainClusters =="1_2"),]$Ann_v4 <- "Unknown_rel2"
NewMeta[which(NewMeta$LouvainClusters =="1_3"),]$Ann_v4 <- "Unknown_rel2"
NewMeta[which(NewMeta$LouvainClusters =="1_4"),]$Ann_v4 <- "Unknown1"
NewMeta[which(NewMeta$LouvainClusters =="2"),]$Ann_v4 <- "PhloemPrecursor"

NewMeta[which(NewMeta$LouvainClusters =="3_1"),]$Ann_v4 <- "Unknown_lowTn5"
NewMeta[which(NewMeta$LouvainClusters =="3_2"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="3_3"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="3_4"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="4_1"),]$Ann_v4 <- "IM-OC"
NewMeta[which(NewMeta$LouvainClusters =="4_2"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="4_3"),]$Ann_v4 <- "Unknown_lowTn5"


NewMeta[which(NewMeta$LouvainClusters =="5"),]$Ann_v4 <- "SPM-base_SM-base"
NewMeta[which(NewMeta$LouvainClusters =="6"),]$Ann_v4 <- ""
NewMeta[which(NewMeta$LouvainClusters =="7"),]$Ann_v4 <- ""
NewMeta[which(NewMeta$LouvainClusters =="8"),]$Ann_v4 <- ""
NewMeta[which(NewMeta$LouvainClusters =="9"),]$Ann_v4 <- ""
NewMeta[which(NewMeta$LouvainClusters =="10"),]$Ann_v4 <- ""
NewMeta[which(NewMeta$LouvainClusters =="11"),]$Ann_v4 <- ""
NewMeta[which(NewMeta$LouvainClusters =="12"),]$Ann_v4 <- ""

levels(factor(NewMeta$Ann_v4))
table(NewMeta$Ann_v4)
colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#0bd43d","#fc53b6",
            "#deadce","#adafde","#5703ff")
length(colorr)

ggplot(NewMeta, aes(x=umap1, y=umap2, color=factor(Ann_v4))) +
  geom_point(size=0.02) +
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))

ggsave("Ref_AnnV4.pdf", width=13, height=10)

write.table(NewMeta, file=paste0("Ref_AnnV4_metadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")

