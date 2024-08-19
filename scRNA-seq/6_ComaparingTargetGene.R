library(Seurat)
library(ggplot2)

obj_WT <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_WTRe1andRe2_UMI1000.rds")
head(obj_WT)
obj_Bif3 <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_Bif3Re1andRe2_UMI1000.rds")

## cluster *6 or cluster 10 - WT CZ
## Cluster 2 and *9 - Bif3 CZ
ClusterTable_WT <- data.frame(obj_WT$seurat_clusters)
ClusterTable_Bif3 <- data.frame(obj_Bif3$seurat_clusters)

CZCells_WT_6 <- rownames(ClusterTable_WT)[which(ClusterTable_WT$obj_WT.seurat_clusters == 6)]
CZCells_WT_10 <- rownames(ClusterTable_WT)[which(ClusterTable_WT$obj_WT.seurat_clusters == 10)]
CZCells_Bif3_2 <- rownames(ClusterTable_Bif3)[which(ClusterTable_Bif3$obj_Bif3.seurat_clusters == 2)]
CZCells_Bif3_9 <- rownames(ClusterTable_Bif3)[which(ClusterTable_Bif3$obj_Bif3.seurat_clusters == 9)]

WT_NomralizedValues <- GetAssayData(object = obj_WT, layer = "data") # Normalized data matrix
head(WT_NomralizedValues[,c(1:10)])
#WT_NomralizedValues["Zm00001eb430980",CZCells_WT_6]
Bif3_NomralizedValues <- GetAssayData(object = obj_Bif3, layer = "data") # Normalized data matrix

###################################
AllPlotgenes_symbol <- c('arftf4', 'arftf30', 'arftf18', 'arftf3', 'arftf20', 'arftf25', 'arftf10', 'arftf36', 'arftf23', 'arftf26', 'knox1')
AllPlotgenes <- c('Zm00001eb433460', 'Zm00001eb066640', 'Zm00001eb067270', 'Zm00001eb142540', 'Zm00001eb224680', 'Zm00001eb232120', 'Zm00001eb243930', 'Zm00001eb292830', 'Zm00001eb363810', 'Zm00001eb370810', 'Zm00001eb001720')
Plotlist <- list()
i <-11
for (i in c(1:length(AllPlotgenes_symbol))) {
  # Extract gene data for WT and Bif3
  GeneID <- AllPlotgenes[[i]]
  GeneName <- AllPlotgenes_symbol[[i]]
  Plot_data <- data.frame(Group="WT_C6", Value=WT_NomralizedValues[GeneID,CZCells_WT_6])
  Plot_data <-rbind(Plot_data, data.frame(Group="WT_C10", Value=WT_NomralizedValues[GeneID,CZCells_WT_10]))
  Plot_data <-rbind(Plot_data, data.frame(Group="Bif3_C2", Value=Bif3_NomralizedValues[GeneID,CZCells_Bif3_2]))
  Plot_data <-rbind(Plot_data, data.frame(Group="Bif3_C9", Value=Bif3_NomralizedValues[GeneID,CZCells_Bif3_9]))
  p <- ggplot(Plot_data, aes(x = Group, y = Value)) +
    geom_boxplot() +
    labs(title = GeneName, x = "Group", y = "Expression Level") +
    theme_minimal()
  Plotlist[[i]] <- p
}

library(cowplot)

final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
output_name <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/6.TargetGene/BoxPlot_11Genes.pdf"
ggsave(output_name, plot = final_plot,
       width = 15, height = 4,
       units = c('in'), limitsize = FALSE,
       dpi = 300)

#### Violet plot ##################
Plotlist <- list()
i <-1
for (i in c(1:length(AllPlotgenes_symbol))) {
  # Extract gene data for WT and Bif3
  GeneID <- AllPlotgenes[[i]]
  GeneName <- AllPlotgenes_symbol[[i]]
  Plot_data <- data.frame(Group="WT_C6", Value=WT_NomralizedValues[GeneID,CZCells_WT_6])
  Plot_data <-rbind(Plot_data, data.frame(Group="Bif3_C9", Value=Bif3_NomralizedValues[GeneID,CZCells_Bif3_9]))
  p <- ggplot(Plot_data, aes(x = Group, y = Value)) +
    geom_boxplot() +
    geom_violin(trim = FALSE) +
    geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
    labs(title = GeneName, x = "Group", y = "Expression Level") +
    theme_minimal()
  Plotlist[[i]] <- p
}

library(cowplot)

final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
output_name <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/6.TargetGene/VioletPlot_11Genes.pdf"
ggsave(output_name, plot = final_plot,
       width = 12, height = 4,
       units = c('in'), limitsize = FALSE,
       dpi = 300)


#### Dot plot (Percentage expressed/Scaled average expression ) or heatmap ######################
Plotlist <- list()
i <-1
for (i in c(1:length(AllPlotgenes_symbol))) {
  # Extract gene data for WT and Bif3
  GeneID <- AllPlotgenes[[i]]
  GeneName <- AllPlotgenes_symbol[[i]]
  print(GeneName)
  print(mean(WT_NomralizedValues[GeneID,CZCells_WT_6]))
  print(mean(Bif3_NomralizedValues[GeneID,CZCells_Bif3_9]))
  print("NonZero")
  print(sum(WT_NomralizedValues[GeneID,CZCells_WT_6]!=0) / length(WT_NomralizedValues[GeneID,CZCells_WT_6]))
  print(sum(Bif3_NomralizedValues[GeneID,CZCells_Bif3_9]!=0) / length(Bif3_NomralizedValues[GeneID,CZCells_Bif3_9]))
  
}

#### Statistical Test ####################################################################
