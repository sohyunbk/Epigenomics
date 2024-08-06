
setwd("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/5.MarkerGene")
Markers <- read.table("/scratch/sb14489/3.scATAC/0.Data/MarkerGene/SpatialMarkerFinal.txt",head=TRUE)
head(Markers)
Markers$geneID[1:10]

MarkerGenePlot <- 
  FeaturePlot(Harmony, features = Markers$geneID[1:10])

ggsave(plot=MarkerGenePlot,"Temp.pdf", width = 10, height = 20)

########################################
library(impute)
library(destiny)

gene_of_interest <- "Zm00001eb067310"

# Assuming 'seurat_obj' is your Seurat object
data_matrix <- GetAssayData(Harmony, slot = "data")

# Create a diffusion map
dm <- DiffusionMap(data_matrix)

# Extract smoothed data (example using first component)
smoothed_data <- dm@eigenvectors[, 1]

# Update Seurat object with smoothed data
Harmony <- SetAssayData(Harmony, slot = "data", new.data = smoothed_data)

# Visualize with FeaturePlot
MarkerGenePlot <- 
   FeaturePlot(Harmony, features = "Zm00001eb067310")
ggsave(plot=MarkerGenePlot,"Temp.pdf", width = 4, height = 3)
