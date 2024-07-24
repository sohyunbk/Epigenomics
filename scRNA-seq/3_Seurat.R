library(stringi)
library(SeuratData)
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)

##############################################################
#############  1) load RNA-seqdata #########################
#############################################################

data_Re1 <- Read10X(data.dir = "/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/2.Mapped/B73re1/outs/filtered_feature_bc_matrix/")
data_Re2 <- Read10X(data.dir = "/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/2.Mapped/B73re2/outs/filtered_feature_bc_matrix/")
data_Re3 <- Read10X(data.dir = "/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/2.Mapped/B73re3/outs/filtered_feature_bc_matrix/")
CellAnnotationFile <- "/scratch/sb14489/4.scRNAseq/0.Data/Dev_Cell_maize_scRNAseq_MetaData_Edited.txt"
Ann <- read.table(CellAnnotationFile,header=TRUE)
head(Ann)
levels(factor(Ann$SampleID))
#"SI-GA-F8" "SI-GA-G3" "SI-GA-G9"
Ann_F8 <- Ann[which(Ann$SampleID=="SI-GA-F8"),]
dim(Ann_F8)
Ann_G3 <- Ann[which(Ann$SampleID=="SI-GA-G3"),]
dim(Ann_G3)
Ann_G9 <- Ann[which(Ann$SampleID=="SI-GA-G9"),]
dim(Ann_G9)
library(dplyr)
table(as.factor(Ann$Celltype))

## MATCHING samples
sum(colnames(data_Re1)%in%Ann_F8$Barcode)
sum(colnames(data_Re1)%in%Ann_G3$Barcode)
sum(colnames(data_Re1)%in%Ann_G9$Barcode)
sum(Ann_F8$Barcode%in%colnames(data_Re1))
## Re1 == F8 : 2250
sum(colnames(data_Re2)%in%Ann_F8$Barcode)
sum(colnames(data_Re2)%in%Ann_G3$Barcode)
sum(colnames(data_Re2)%in%Ann_G9$Barcode)
sum(Ann_G3$Barcode%in%colnames(data_Re2))

## Re2 == G3 : 6192

sum(colnames(data_Re3)%in%Ann_F8$Barcode)
sum(colnames(data_Re3)%in%Ann_G3$Barcode)
sum(colnames(data_Re3)%in%Ann_G9$Barcode)
sum(Ann_G9$Barcode%in%colnames(data_Re3))
## Re3 == G9 : 4076
2250+6192+4076
dim(Ann)

########### Create Seurat object & Summarize the matrix 
data_Re1_Filtered <- data_Re1[,colnames(data_Re1)%in%Ann_F8$Barcode]
data_Re2_Filtered <- data_Re2[,colnames(data_Re2)%in%Ann_G3$Barcode]
data_Re3_Filtered <- data_Re3[,colnames(data_Re3)%in%Ann_G9$Barcode]
Ann_F8$Replicates <- "Re1"
Ann_G3$Replicates <- "Re2"
Ann_G9$Replicates <- "Re3"
MetaData <- rbind(Ann_F8,Ann_G3,Ann_G9)
head(MetaData)
MetaData$newBarcode <- paste(tolower(MetaData$Replicates),MetaData$Barcode,
                             MetaData$Replicates,sep="_")

tail(MetaData)
write.table(MetaData,"Metadata_Edit.txt",col.names=TRUE)

## Change the column name by adding _Re1
colnames(data_Re1_Filtered) <- paste0(colnames(data_Re1_Filtered),"_Re1")
colnames(data_Re2_Filtered) <- paste0(colnames(data_Re2_Filtered),"_Re2")
colnames(data_Re3_Filtered) <- paste0(colnames(data_Re3_Filtered),"_Re3")


Obj_re1 <- CreateSeuratObject(counts = data_Re1_Filtered, project = "Re1", min.cells = 3, min.features = 200)
#Obj_re1 <- subset(Obj_re1, cells=Ann_F8$Barcode)
str(Obj_re1)
Obj_re2 <- CreateSeuratObject(counts = data_Re2_Filtered, project = "Re2", min.cells = 3, min.features = 200)
str(Obj_re2)
Obj_re3 <- CreateSeuratObject(counts = data_Re3_Filtered, project = "Re3", min.cells = 3, min.features = 200)
str(Obj_re3)

TotalData <- merge(x = Obj_re1, y = c(Obj_re2, Obj_re3),
                 add.cell.ids = c("re1", "re2", "re3"), project = "Maize_ear")
str(TotalData)
25996*0.97
rna <- NormalizeData(TotalData ,normalization.method = "LogNormalize")
rna <- FindVariableFeatures(rna, nfeatures = 25000)
rna <- ScaleData(rna)
str(rna)
length(levels(factor(rownames(rna))))
rna <- RunPCA(rna)
?RunUMAP
rna <- RunUMAP(rna, dims = 1:30)

setwd("/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/4.PublishedCells")
saveRDS(rna,"Allprocessed_Seurat.rds")
saveRDS(TotalData,"AllJustCombined_Seurat.rds")

## Check UMAP
setwd("/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/4.PublishedCells")

str(rna@reductions$umap)
head(rna@reductions$umap@cell.embeddings)
UMAP <- rna@reductions$umap@cell.embeddings
head(MetaData)
dim(UMAP)
dim(MetaData)
MetaData <- MetaData[MetaData$newBarcode%in%rownames(UMAP),]
head(UMAP[MetaData$newBarcode,])
MetaData$umap1 <- UMAP[MetaData$newBarcode,1]
MetaData$umap2 <- UMAP[MetaData$newBarcode,2]

head(MetaData)
levels(factor(MetaData$Celltype))
length(which(MetaData$Celltype=="Cells_unique_to_only_one_replicate"))
length(which(MetaData$Celltype=="Phloem"))
dim(MetaData)

colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#de62b9","#0bd43d","#cf1d35","#deadce","#adafde")

ggplot(MetaData, aes(x=umap1, y=umap2, color=factor(Celltype))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("All_withCelltype_PC30.pdf", width=13, height=10)	

MetaData_o <- MetaData[MetaData$Celltype!="Cells_unique_to_only_one_replicate",]
head(MetaData_o)
ggplot(MetaData_o, aes(x=umap1, y=umap2, color=factor(Celltype))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("OnlyValidCells_withCelltype_PC30.pdf", width=13, height=10)	


ggplot(MetaData_o, aes(x=umap1, y=umap2, color=factor(Replicates))) +
  geom_point(size=0.02) + 
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("OnlyValidCells_Replicates_PC30.pdf", width=13, height=10)	


ggplot(MetaData, aes(x=umap1, y=umap2, color=factor(Replicates))) +
  geom_point(size=0.02) + 
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("All_Replicates_PC30.pdf", width=13, height=10)


str(rna)
str(rna@assays)
head(rna@assays$RNA@counts)

##############################################################
#############  2) load ATAC-seqdata #########################
#############################################################
ATAC <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.afterHarmony.processed.rds")
str(ATAC)
obj_ATAC <- ATAC$counts
head(ATAC$counts)[,c(1:5)]
obj_ATAC <- CreateSeuratObject(obj_ATAC, assay = "ATAC")
str(obj_ATAC)
#obj_ATAC



##############################################################
#############  3) RNA - ATAC ################################
#############################################################
gene.activities <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/GA_A619.txt",header=TRUE)
head(gene.activities)
obj_ATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)












#gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))
# add gene activities as a new assay
obj_ATAC[["ACTIVITY"]] <- CreateAssayObject(counts = ATAC$counts)

# normalize gene activities
DefaultAssay(obj_ATAC) <- "ACTIVITY"
obj_ATAC <- NormalizeData(obj_ATAC)
obj_ATAC <- ScaleData(obj_ATAC, features = rownames(obj_ATAC))

transfer.anchors <- FindTransferAnchors(reference = rna, query = obj_ATAC,
                                        features = VariableFeatures(object = rna),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
                                     weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)




# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
