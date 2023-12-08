## Make metafile & plot for estimated annotation
library(ggplot2)
library(stringr)

## Make new metafile bycombining all the meta
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3")
meta <- "bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt"
loaded_meta_data <- read.table(meta)
CellOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order.txt")
Cluster1 <- read.table("bif3_Cluster1_Recluster_Sub_res1_knear100_Partmetadata.txt")
head(Cluster1)


head(loaded_meta_data)
#levels(factor(loaded_meta_data$LouvainClusters))
loaded_meta_data[rownames(Cluster1),]$LouvainClusters <- Cluster1$LouvainClusters
levels(factor(loaded_meta_data$LouvainClusters))

NewMeta <- loaded_meta_data
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="1"),]

dim(loaded_meta_data)
dim(NewMeta)

levels(factor(NewMeta$LouvainClusters))
NewMeta$Ann_v4 <- "Ann_v4"
levels(factor(NewMeta$Ann_v4))
## New assignment ##
NewMeta[which(NewMeta$LouvainClusters =="1_1"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="1_2"),]$Ann_v4 <- "SPM-base_SM-base"
NewMeta[which(NewMeta$LouvainClusters =="1_3"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="2"),]$Ann_v4 <- "Unknown1"
NewMeta[which(NewMeta$LouvainClusters =="3"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="4"),]$Ann_v4 <- "XylemParenchyma_PithParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="5"),]$Ann_v4 <- "ProcambialMeristem_ProtoXylem_MetaXylem"
NewMeta[which(NewMeta$LouvainClusters =="6"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="7"),]$Ann_v4 <- "Unknown2"
NewMeta[which(NewMeta$LouvainClusters =="8"),]$Ann_v4 <- "L1atFloralMeristem"
NewMeta[which(NewMeta$LouvainClusters =="9"),]$Ann_v4 <- "L1"
NewMeta[which(NewMeta$LouvainClusters =="10"),]$Ann_v4 <- "Unknown_lowFRiP"
NewMeta[which(NewMeta$LouvainClusters =="11"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="12"),]$Ann_v4 <- "G2_M"
NewMeta[which(NewMeta$LouvainClusters =="13"),]$Ann_v4 <- "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="14"),]$Ann_v4 <- "IM-OC"
######################

levels(factor(NewMeta$Ann_v4))
table(NewMeta$Ann_v4)
colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#0bd43d","#fc53b6",
            "#deadce","#adafde","#5703ff")
length(colorr)
NewMeta$Ann_v4 <- factor(NewMeta$Ann_v4, levels = CellOrder)

ggplot(NewMeta, aes(x=umap1, y=umap2, color=factor(Ann_v4))) +
  geom_point(size=0.02) +
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme(text=element_text(size=20))

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3")
ggsave("Bif3_AnnV4.pdf", width=14, height=8)

write.table(NewMeta, file=paste0("Bif3_AnnV4_metadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")
