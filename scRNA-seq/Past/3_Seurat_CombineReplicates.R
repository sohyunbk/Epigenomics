# conda activate r_env_copy
library(stringi)
library(SeuratData)
#setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
#install.packages("Signac")
# install the dataset and load requirements
#InstallData("pbmcMultiome")
library(Seurat)
#library(Signac)

#library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(stringr)

#########################
# Load the dataset #
#########################

setwd("/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/3.Seurat")

#################################
## load Re1
Re1 <- Read10X(data.dir = "/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/2.Mapped/B73re1/outs/filtered_feature_bc_matrix/")
str(Re1)
# Gene	39756

Re1[c(1:20),c(1:20)]
Obj_re1 <- CreateSeuratObject(counts = Re1, project = "Re1", min.cells = 3, min.features = 200)
str(Obj_re1)
head(Obj_re1@meta.data)
# Re1: 3304
# 
Obj_re1[["percent.mt"]] <- PercentageFeatureSet(Obj_re1, pattern = "GRM") ## Annotation of gtf
table(Obj_re1[["percent.mt"]])
head(Obj_re1@meta.data, 5)

pdf(file="Re1_VioletPlot_nFeature_nCount_MtPer.pdf")
VlnPlot(Obj_re1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Obj_re1 <- subset(Obj_re1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 2)
str(Obj_re1)
#?PercentageFeatureSet
## load Re2
Re2 <- Read10X(data.dir = "/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/2.Mapped/B73re2/outs/filtered_feature_bc_matrix/")
Re2[c(1:20),c(1:20)]
str(Re2)

Obj_re2 <- CreateSeuratObject(counts = Re2, project = "Re2", min.cells = 3, min.features = 200)
str(Obj_re2)
Obj_re2[["percent.mt"]] <- PercentageFeatureSet(Obj_re2, pattern = "GRM") ## Annotation of gtf
table(Obj_re2[["percent.mt"]])
head(Obj_re2@meta.data, 5)


pdf(file="Re2_VioletPlot_nFeature_nCount_MtPer.pdf")
VlnPlot(Obj_re2, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf(file="Re2_ScatterPlots_forQC.pdf",width=20)
plot1 <- FeatureScatter(Obj_re2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Obj_re2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


Obj_re2 <- subset(Obj_re2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 2)
str(Obj_re2)

## Road Re3
Re3 <- Read10X(data.dir = "/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/2.Mapped/B73re3/outs/filtered_feature_bc_matrix/")
str(Re3)

Obj_re3 <- CreateSeuratObject(counts = Re3, project = "Re3", min.cells = 3, min.features = 200)
str(Obj_re3)

Obj_re3[["percent.mt"]] <- PercentageFeatureSet(Obj_re3, pattern = "GRM") ## Annotation of gtf
table(Obj_re3[["percent.mt"]])
head(Obj_re3@meta.data, 5)

pdf(file="Re3_VioletPlot_nFeature_nCount_MtPer.pdf")
VlnPlot(Obj_re3, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

pdf(file="Re3_ScatterPlots_forQC.pdf",width=20)
plot1 <- FeatureScatter(Obj_re3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Obj_re3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


Obj_re3 <- subset(Obj_re3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 2)
str(Obj_re2)



#########################
# Merge dataset
#########################
# merge two objects
# to merge more than two objects, pass one to x and a list of objects to y
obj_com <- merge(x = Obj_re1, y = c(Obj_re2, Obj_re3),add.cell.ids = c("re1", "re2", "re3"), project = "Maize_ear")
#https://satijalab.org/seurat/articles/merge_vignette.html
#Re1	3304
#Re2	9925
#Re3	5880
#Sum	19109

str(obj_com)
head(obj_com@meta.data)
tail(obj_com@meta.data)

##Normalizing the data
obj_com <- NormalizeData(obj_com, normalization.method = "LogNormalize", scale.factor = 10000)
#obj_com <- NormalizeData(obj_com)

## Identification of highly variable features (feature selection)
obj_com <- FindVariableFeatures(obj_com, selection.method = "vst", nfeatures = 8000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj_com), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj_com)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file="HighlyVariableFeatures.pdf",width=20)
CombinePlots(plots = list(plot1, plot2))
dev.off()

## Scaling the data
## Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function:

all.genes <- rownames(obj_com)
obj_com <- ScaleData(obj_com, features = all.genes)
head(obj_com)
str(obj_com)
## Perform linear dimensional reduction
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
obj_PCA <- RunPCA(obj_com, features = VariableFeatures(object = obj_com))
# Examine and visualize PCA results a few different ways
print(obj_PCA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(obj_PCA, dims = 1:2, reduction = "pca")

pdf(file="PCA.pdf",width=10)
DimPlot(obj_PCA, reduction = "pca")
dev.off()

pdf(file="Heatmap_PC1.pdf",width=15)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf(file="Heatmap_PC1_15.pdf",width=15)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

?PercentageFeatureSet
data.filt_v <- FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = 8000)








all.genes <- rownames(data.filt_v )
data.filt_v  <- ScaleData(data.filt_v , features = all.genes)
str(data.filt_v )
head(data.filt_v [["RNA"]]@scale.data)


data.filt_v  <- RunPCA(data.filt_v , features = VariableFeatures(object = data.filt_v ))
print(data.filt_v[["pca"]], dims = 1:5, nfeatures = 5)


PCA <- FindNeighbors(data.filt_v, k.param=,dims = 1:10)
PCA <- FindClusters(PCA, resolution = 0.5)


UMAP <- RunUMAP(PCA, dims = 1:10)

str(UMAP)

setwd("/scratch/sb14489/4.scRNAseq/")
tiff(file=paste("AllReplicates_UMAP.tiff",sep=""),type="cairo")
DimPlot(UMAP, reduction = "umap")
dev.off()

# How do I create a UMAP plot where cells are colored by replicate?  First, store the current
# identities in a new column of meta.data called CellType


# Next, switch the identity class of all cells to reflect replicate ID
Idents(pbmc) <- c(rep("Re1",ncol(Re1)),rep("Re2",ncol(Re2)),rep("Re3",ncol(Re3)))
tiff(file=paste("AllReplicates_UMAP.tiff",sep=""),type="cairo")
DimPlot(pbmc, reduction = "umap")
dev.off()
####################################
