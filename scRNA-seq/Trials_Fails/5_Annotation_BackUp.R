Harmony <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_WTRe1andRe2.rds")
Name <- "WTRe1andRe2"

setwd("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/5.MarkerGene")
Markers <- read.table("/scratch/sb14489/3.scATAC/0.Data/MarkerGene/SeletedMarkergeneForDotPlot_RemoveUnknown.txt",head=TRUE)
head(Markers)
dim(Markers)

MarkerGenePlot <- 
  FeaturePlot(Harmony, features = Markers$geneID)

ggsave(plot=MarkerGenePlot,paste0(Name,"_28MarkerGene.pdf"), width = 10, height = 22)


############################





Harmony <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_WTRe1andRe2.rds")
Harmony <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_Bif3Re1andRe2.rds")

MarkerGenePlot <- FeaturePlot(Harmony, features = "Zm00001eb067310")
ggsave(plot=MarkerGenePlot,"Bif3Re1andRe2_ZmWUS1.pdf", width = 10, height = 10)

normalized_data <- GetAssayData(obj, layer = "data")
raw_data <- GetAssayData(obj, slot = "counts")
raw_data2 <- GetAssayData(obj_filtered, slot = "counts")

sum(GetAssayData(Harmony, layer = "counts")["Zm00001eb999999",])
sum(GetAssayData(Harmony, layer = "counts")["Zm00001eb067310",])





########################################
library(tradeSeq)
library(Seurat)

obj <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/obj_afterDoubletWT_Re2.rds")
obj <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/obj_afterDoubletWT_Re1.rds")

obj <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/obj_afterDoubletBif3_Re1.rds")
obj <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat/obj_afterDoubletBif3_Re2.rds")

normalized_data <- GetAssayData(obj, layer = "data")
raw_data <- GetAssayData(obj, slot = "counts")
raw_data2 <- GetAssayData(obj_filtered, slot = "counts")

sum(normalized_data["Zm00001eb999999",])

sum(raw_data["Zm00001eb067310",])
sum(raw_data2["Zm00001eb067310",])

sum(raw_data["Zm00001eb999999",])

obj_filtered

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

