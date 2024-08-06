library(harmony)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(gridExtra)
library(grid)  # Load the grid package
library("RColorBrewer")
library("optparse")

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name", metavar="character"),
  make_option(c("--rds1"), type="character",
              help="rds1", metavar="character"),
  make_option(c("--rds2"), type="character",
              help="rds2", metavar="character")
 
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd(opt$WD)
#setwd("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony")

DataName <- opt$Name
Re1 <- opt$rds1
Re2 <- opt$rds2

Re1and2 <- merge(Re1,y=Re2,add.cell.ids=c("Re1","Re2"))
Re1and2 <- JoinLayers(Re1and2)
Re1and2@meta.data$Replicates <- sapply(strsplit(rownames(Re1and2@meta.data), "_"), function(x) x[1])

obj_filtered <- NormalizeData(Re1and2)
obj_filtered <- FindVariableFeatures(obj_filtered, selection.method = "vst", nfeatures = 4000)
obj_filtered <- ScaleData(obj_filtered)
obj_filtered <- RunPCA(obj_filtered)

Harmony <- RunHarmony(obj_filtered, group.by.vars = "Replicates")
Harmony <- RunUMAP(Harmony, reduction = "harmony", dims = 1:10)
Harmony <- FindNeighbors(Harmony, reduction = "harmony", dims = 1:10)
Harmony <- FindClusters(Harmony, resolution = 1)
saveRDS(Harmony, file = paste0("obj_afterHarmony_",DataName,".rds"))

# check embeddings
pca_embeddings <- Embeddings(Harmony, "pca")
harmony_embeddings <- Embeddings(Harmony, "harmony")
head(pca_embeddings)
head(harmony_embeddings)


#Harmony_umap_plot <- DimPlot(Harmony, reduction = "umap")

NoHarmony <- FindNeighbors(obj_filtered, dims = 1:10)
NoHarmony <- FindClusters(NoHarmony, resolution = 1)
NoHarmony <- RunUMAP(NoHarmony, dims = 1:10)
NoHarmony_umap_plot <- DimPlot(NoHarmony, reduction = "umap")
###############################################################
Meta <- data.frame(Harmony@meta.data)
UMAPTable <- data.frame(Harmony@reductions$umap@cell.embeddings)
ClusterTable <- data.frame(Harmony$seurat_clusters)

Meta$cell_id <- rownames(Meta)
UMAPTable$cell_id <- rownames(UMAPTable)
ClusterTable$cell_id <- rownames(ClusterTable)

UMAPTable <- merge(Meta, UMAPTable, by = "cell_id")
UMAPTable <- merge(UMAPTable, ClusterTable, by = "cell_id")

dim(UMAPTable)
head(UMAPTable)


### 1) Cluster
colorr <- c("#4F96C4","#84f5d9","#0bd43d","#d62744","#FDA33F","#060878","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#adafde","#DE9A89","#5703ff",
            "#deadce","#fc53b6")

All <- ggplot(UMAPTable, aes(x=umap_1, y=umap_2, color=factor(Harmony.seurat_clusters))) +
  geom_point(size=0.05) +
  scale_color_manual(values=colorr)+
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=12)))+
  labs(title = paste0("Re1+R2 \n CellNumber: ",nrow(UMAPTable)),
       x = "UMAP1",
       y = "UMAP2")+
  theme(legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

UMAPReplicate <- ggplot(UMAPTable, aes(x=umap_1, y=umap_2, color=factor(Replicates))) +
  geom_point(size=0.05) +
  scale_color_manual(values=c("red","blue"))+
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=12)))+
  labs(title = paste0("Re1+R2 \n CellNumber: ",nrow(UMAPTable)),
       x = "UMAP1",
       y = "UMAP2")+
  theme(legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))


UMAPUMICount <- ggplot(UMAPTable, aes(x = umap_1, y = umap_2, color = log10_nCount_RNA)) +
  geom_point(size = 0.05) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("blue", "grey", "red"))(100),
    limits = c(min(UMAPTable$log10_nCount_RNA), max(UMAPTable$log10_nCount_RNA))
  ) +
  theme_minimal() +
  labs(title = "UMI Counts",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "log10(nCount_RNA)")

UMAP_Org <- ggplot(UMAPTable, aes(x = umap_1, y = umap_2, color = percent.organelle)) +
  geom_point(size = 0.05) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("blue", "grey", "red"))(100),
    limits = c(min(UMAPTable$percent.organelle), max(UMAPTable$percent.organelle))
  ) +
  theme_minimal() +
  labs(title = "percent.organelle",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "Org %")

#ggsave(plot=NoHarmony_umap_plot,"Temp.pdf", width = 10, height = 7.5)
combined_plot <- grid.arrange(
  NoHarmony_umap_plot,All, UMAPReplicate, UMAPUMICount,UMAP_Org,
  ncol = 5,
  widths = c(3,4, 3,3,3)
)

# Save the combined plot
ggsave(filename = paste0("Harmony_", DataName, ".pdf"),
       plot = combined_plot, width = 30, height = 7)


