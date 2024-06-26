## Make metafile & plot for estimated annotation
library(ggplot2)
library(stringr)
source("/home/sb14489/Epigenomics/scATAC-seq/0_Function/DrawFigures_QC_Annotation_forUMAP.R")

## Make new metafile bycombining all the meta
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/rel2/")
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt"
loaded_meta_data <- read.table(meta)
PreAnn <- loaded_meta_data
head(loaded_meta_data)

CellOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forrel2.txt")

#For rel2, I did sub clustering for cluster 1,3,4
SubClusterDir <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/4.Subclustering/rel2/"
Cluster1 <- read.table(paste0(SubClusterDir,"Cluster1_Sub_res1_knear100_Partmetadata.txt"))
Cluster1$LouvainClusters <- paste0("1_",Cluster1$LouvainClusters)
Cluster2 <- read.table(paste0(SubClusterDir,"Cluster2_Sub_res2_knear100_Partmetadata.txt"))
head(Cluster2)
Cluster2$LouvainClusters <- paste0("2_",Cluster2$LouvainClusters)
Cluster3 <- read.table(paste0(SubClusterDir,"Cluster3_Sub_res1_knear100_Partmetadata.txt"))
Cluster3$LouvainClusters <- paste0("3_",Cluster3$LouvainClusters)
Cluster4 <- read.table(paste0(SubClusterDir,"Cluster4_Sub_res1_knear100_Partmetadata.txt"))
Cluster4$LouvainClusters <- paste0("4_",Cluster4$LouvainClusters)
head(Cluster1)
#levels(factor(loaded_meta_data$LouvainClusters))
loaded_meta_data[rownames(Cluster1),]$LouvainClusters <- Cluster1$LouvainClusters
loaded_meta_data[rownames(Cluster2),]$LouvainClusters <- Cluster2$LouvainClusters
loaded_meta_data[rownames(Cluster3),]$LouvainClusters <- Cluster3$LouvainClusters
loaded_meta_data[rownames(Cluster4),]$LouvainClusters <- Cluster4$LouvainClusters

levels(factor(loaded_meta_data$LouvainClusters))

NewMeta <- loaded_meta_data
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="1"),]
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="3"),]
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="4"),]
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="2"),]

dim(loaded_meta_data)
dim(NewMeta)
levels(factor(NewMeta$LouvainClusters))

NewMeta$Ann_v4 <- "Ann_v4"
levels(factor(NewMeta$Ann_v4))

NewMeta[which(NewMeta$LouvainClusters =="1_1"),]$Ann_v4 <- "Unknown1"
NewMeta[which(NewMeta$LouvainClusters =="1_2"),]$Ann_v4 <- "Unknown_rel2"
NewMeta[which(NewMeta$LouvainClusters =="1_3"),]$Ann_v4 <- "Unknown_rel2"
NewMeta[which(NewMeta$LouvainClusters =="1_4"),]$Ann_v4 <- "Unknown1"

NewMeta[which(NewMeta$LouvainClusters =="2_1"),]$Ann_v4 <- "G2_M"
NewMeta[which(NewMeta$LouvainClusters =="2_2"),]$Ann_v4 <- "ProcambialMeristem_ProtoXylem_MetaXylem"
NewMeta[which(NewMeta$LouvainClusters =="2_3"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="2_4"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="2_5"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="2_6"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="2_7"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="2_8"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="2_9"),]$Ann_v4 <- "PhloemPrecursor"

NewMeta[which(NewMeta$LouvainClusters =="3_1"),]$Ann_v4 <- "Unknown_lowFRiP"
NewMeta[which(NewMeta$LouvainClusters =="3_2"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="3_3"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="3_4"),]$Ann_v4 <- "Unknown_Sclerenchyma"

NewMeta[which(NewMeta$LouvainClusters =="4_1"),]$Ann_v4 <- "IM-OC"
NewMeta[which(NewMeta$LouvainClusters =="4_2"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="4_3"),]$Ann_v4 <- "Unknown_lowFRiP"


NewMeta[which(NewMeta$LouvainClusters =="5"),]$Ann_v4 <- "SPM-base_SM-base"
NewMeta[which(NewMeta$LouvainClusters =="6"),]$Ann_v4 <- "L1"
NewMeta[which(NewMeta$LouvainClusters =="7"),]$Ann_v4 <- "L1atFloralMeristem"
NewMeta[which(NewMeta$LouvainClusters =="8"),]$Ann_v4 <- "XylemParenchyma_PithParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="9"),]$Ann_v4 <- "XylemParenchyma_PithParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="10"),]$Ann_v4 <- "G2_M"
NewMeta[which(NewMeta$LouvainClusters =="11"),]$Ann_v4 <- "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="12"),]$Ann_v4 <- "ProcambialMeristem_ProtoXylem_MetaXylem"

write.table(NewMeta, file=paste0("rel2_AnnV4.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")
DrawUMAP_Ann_QC(PreAnn,NewMeta, "Ann_v4", CellOrder, "rel2_Re1", "rel2_Re2","rel2_AnnV4")

