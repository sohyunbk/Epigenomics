library(stringi)
library(Seurat)
#library(Signac)
library(ggplot2)
library(cowplot)
library(tidyr)
library(gridExtra)
library(grid)  # Load the grid package
library(DoubletFinder)
library("optparse")

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name", metavar="character"),
  make_option(c("--UMICut"), type="integer",
              help="UMICut", metavar="integer"),
  make_option(c("--GeneCut"), type="character",
              help="GeneCut", metavar="character"),
  make_option(c("--InputDir"), type="character",
              help="InputDir", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd(opt$WD)

###############################################################################
### 1) Load raw counts from cellRanger
###############################################################################
#DataName <- "WT_Re1_Gene1000_UMI1000"
DataName <- opt$Name

#InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-wt-1/outs/raw_feature_bc_matrix/"
InputFile <- opt$InputDir
  
umi_count_threshold <- opt$UMICut
gene_count_threshold <- opt$GeneCut

#DataName <- "WT_Re2"
#InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-wt-2/outs/raw_feature_bc_matrix/"

#DataName <- "Bif3_Re1"
#InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-bif3-1/outs/raw_feature_bc_matrix/"

#DataName <- "Bif3_Re2"
#InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-bif3-2/outs/raw_feature_bc_matrix/"

data <- Read10X(data.dir = InputFile)
obj <- CreateSeuratObject(counts = data, project = DataName, min.cells = 0, min.features = 0)
obj[["percent.organelle"]] <- PercentageFeatureSet(obj, pattern = "^gene:GRMZM")

obj$log10_nFeature_RNA <- log10(obj$nFeature_RNA + 0.00001)  # Add 1 to avoid log(0)
obj$log10_nCount_RNA <- log10(obj$nCount_RNA + 0.00001)  # Add 1 to avoid log(0)
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

################################################################################
### 2) Knee plot
################################################################################
min.cells <- 1000
max.cells <- 16500
umi <- obj@meta.data[order(obj@meta.data$nCount_RNA, decreasing=T),]
head(umi)
tail(umi)
umi$rank <- 1:nrow(umi)
# The smoothing parameter spar controls the trade-off between smoothness and fidelity to the data. A smaller value of spar makes the fit more flexible (less smooth), while a larger value makes it smoother.
df <- data.frame(rank=log10(umi$rank), depth=log10(umi$nCount_RNA+1))
head(df)
fit <- smooth.spline(df$rank, df$depth, spar=0.1)
X <- data.frame(t = seq(min(df$rank), max(df$rank), length = nrow(df)))
Y <- predict(fit, newdata=X, deriv=1)
str(Y)
xvals <- Y$x
yvals <- Y$y
knee <- xvals[which.min(yvals[min.cells:max.cells])]
cells <- which.min(yvals[min.cells:max.cells])
UMICutoff <- (10^(min(df[1:cells,]$depth))) - 1

obj_AfterKnee <- subset(obj, subset = nCount_RNA > UMICutoff)
ncol(obj)
ncol(obj_AfterKnee)

KneePlot <- ggplot(umi, aes(x = log10(1:nrow(umi)), y = log10(nCount_RNA))) +
  geom_line(size=3) +
  labs(x = "Barcode rank (log10) ", y = "UMI counts (log10) ") +
  ggtitle(paste0(DataName,"\nCell # after filter:",ncol(obj_AfterKnee),
                 " | Cutoff UMI #:",UMICutoff)) +
  geom_hline(yintercept=log10(UMICutoff), linetype="dashed", color = "red") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

#ggsave(filename = paste0("Kneeplot_", DataName, ".pdf"), plot = KneePlot, width = 6, height = 5.5)

### 3-1) Gene count - feature count filter
nfold <- 1
median_gene_count <- median(obj_AfterKnee$nFeature_RNA)
mad_gene_count <- mad(obj_AfterKnee$nFeature_RNA)

#gene_count_threshold <- 1000
if (gene_count_threshold == 0) {
  gene_count_threshold <- median_gene_count - (nfold*mad_gene_count)
}


## 3-2) UMI cutoff
if (umi_count_threshold == 0) {
  umm_count_threshold <- quantile(obj_AfterKnee$nCount_RNA, 0.10)
}

#umi_count_threshold <- quantile(obj_AfterKnee$nCount_RNA, 0.10)
#umi_count_threshold <- 1000
#gene_count_threshold <- 1000

## 3-3) Org cutoff 
nOrg <- 5

## 3-4) Genecount / UMI
obj_filtered <- subset(obj_AfterKnee, subset = nCount_RNA > umi_count_threshold &
                         nFeature_RNA > gene_count_threshold &
                         percent.organelle < nOrg )

head(obj_filtered )
obj_filtered

### violet Plot
p1 <- ggplot(obj_AfterKnee@meta.data, aes(x=orig.ident, y=log10_nFeature_RNA)) + 
  geom_violin()+
  geom_hline(yintercept = log10(gene_count_threshold), linetype = "dashed", color = "red") + 
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Cutoff: ", gene_count_threshold), hjust = 1.1, vjust = 2, size = 5, color = "black")+
  labs(x = DataName, y = "Log10 Gene Count")

p2 <- ggplot(obj_AfterKnee@meta.data, aes(x=orig.ident, y=log10_nCount_RNA)) + 
  geom_violin()+
  geom_hline(yintercept = log10(umi_count_threshold), linetype = "dashed", color = "red") + 
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Cutoff: ", umi_count_threshold), hjust = 1.1, vjust = 2, size = 5, color = "black")+
  labs(x = DataName, y = "Log10 UMI Count")

p3 <- ggplot(obj_AfterKnee@meta.data, aes(x=orig.ident, y=percent.organelle)) + 
  geom_violin()+
  geom_hline(yintercept = nOrg, linetype = "dashed", color = "red") + 
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Cutoff: ", nOrg), hjust = 1.1, vjust = 2, size = 5, color = "black")+
  
  labs(x = DataName, y = "Organelle Percentage")

p4 <- obj_AfterKnee@meta.data %>% 
  ggplot(aes(x=log10GenesPerUMI)) + 
  theme_classic() +
  geom_density(fill = "blue", alpha = 0.5) +
  ylab("Cell density") 

###################################
# Plotting  the QC before doublet.
###################################

combined_plot <- grid.arrange(
  KneePlot, p1, p2, p3,p4,
  ncol = 5,
  widths = c(4, 2, 2, 2,4), 
  top = textGrob(paste0("Cell Number # after GeneCount, UMI Count, Organelle Ratio Filtering: ", 
                        ncol(obj_filtered)),
                 gp = gpar(fontsize = 13))
)


# Save the combined plot
ggsave(filename = paste0("QCPlotKneeViolet_", DataName, ".pdf"),
       plot = combined_plot, width = 12, height = 8)


####################################
# 4) Doublet 
####################################

obj_filtered <- NormalizeData(obj_filtered)
obj_filtered <- FindVariableFeatures(obj_filtered, selection.method = "vst", nfeatures = 4000)
obj_filtered <- ScaleData(obj_filtered)
obj_filtered <- RunPCA(obj_filtered)

obj_filtered <- FindNeighbors(obj_filtered, dims = 1:10)
obj_filtered <- FindClusters(obj_filtered, resolution = 1)

obj_filtered <- RunUMAP(obj_filtered, dims = 1:10)
# Assuming 'obj' is your Seurat object and you have run RunUMAP on it

# Display the first few rows of the UMAP coordinates directly from the object
head(obj_filtered@reductions$umap@cell.embeddings)

umap_plot <- DimPlot(obj_filtered, reduction = "umap")
#ggsave("UMAPTest_res1.pdf", plot = umap_plot, width = 8, height = 6)

nCell = ncol(obj_filtered) #calculate cell number in 1k unit
nExpRate = round(nCell/1000) * 0.008
nExp <- round(nCell * nExpRate)         # expected number of doublets  (0.8% double rate, 8 doublets in 1k cells called)
message('...Cell number: ', nCell, ". Expected doublets: ", nExp, "(", nExpRate, ")")
obj_filtered_Doublet <- doubletFinder(obj_filtered, pN = 0.25, pK = 0.09,
                              nExp = nExp, PCs = 1:10, sct=T)   # find doublet It takes long like 20 min

colnames(obj_filtered_Doublet@meta.data)[ncol(obj_filtered_Doublet@meta.data)] <- "DoubletClass"
colnames(obj_filtered_Doublet@meta.data)[ncol(obj_filtered_Doublet@meta.data)-1] <- "DoubletValue"

Meta <- data.frame(obj_filtered_Doublet@meta.data)
UMAPTable <- data.frame(obj_filtered_Doublet@reductions$umap@cell.embeddings)

Meta$cell_id <- rownames(Meta)
UMAPTable$cell_id <- rownames(UMAPTable)

UMAPTable <- merge(Meta, UMAPTable, by = "cell_id")
head(UMAPTable)
##Plot
library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAPTable$DoubletValue), max(UMAPTable$DoubletValue)))

obj_final <- subset(obj_filtered_Doublet, subset = DoubletClass == "Singlet")

#ncol(obj_final)
DoubletUMAPPlot <- ggplot(UMAPTable, aes(x=umap_1, y=umap_2,color=DoubletValue)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc+
  annotate("text", x = Inf, y = Inf,
           label = paste0("Cell # filtering doublet: \n",
                          ncol(obj_final))
                          , hjust = 1.1, vjust = 2, size = 4, color = "black")



##################################
# 5) Other QC plot on UMAP 
##################################
Meta <- data.frame(obj_final@meta.data)
UMAPTable <- data.frame(obj_final@reductions$umap@cell.embeddings)

Meta$cell_id <- rownames(Meta)
UMAPTable$cell_id <- rownames(UMAPTable)

UMAPTable <- merge(Meta, UMAPTable, by = "cell_id")
dim(UMAPTable)
head(UMAPTable)


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAPTable$log10_nFeature_RNA), max(UMAPTable$log10_nFeature_RNA)))
GeneUMAPPlot <- ggplot(UMAPTable, aes(x=umap_1, y=umap_2,color=log10_nFeature_RNA)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc

sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAPTable$log10_nCount_RNA), max(UMAPTable$log10_nCount_RNA)))
UMI_UMAPPlot <- ggplot(UMAPTable, aes(x=umap_1, y=umap_2,color=log10_nCount_RNA)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc
####
head(UMAPTable)
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAPTable$percent.organelle), max(UMAPTable$percent.organelle)))
Org_UMAPPlot <- ggplot(UMAPTable, aes(x=umap_1, y=umap_2,color=percent.organelle)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc
###############################
# Plotting  the QC and save object
###################################

saveRDS(obj_final, file = paste0("obj_afterDoublet",DataName,".rds"))

combined_plot <- grid.arrange(
  umap_plot,DoubletUMAPPlot, GeneUMAPPlot, UMI_UMAPPlot,Org_UMAPPlot,
  ncol = 5,
  widths = c(2,2, 2, 2, 2)
)

# Save the combined plot
ggsave(filename = paste0("QCPlotUMAPDoublet_", DataName, ".pdf"),
       plot = combined_plot, width = 17, height = 7)
