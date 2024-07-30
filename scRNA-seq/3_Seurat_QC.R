library(stringi)
library(SeuratData)
library(Seurat)
#library(Signac)
library(ggplot2)
library(cowplot)
library(tidyr)
library(gridExtra)
library(grid)  # Load the grid package
library(DoubletFinder)

setwd("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat")

###############################################################################
### 1) Load raw counts from cellRanger
###############################################################################
DataName <- "WT_Re1"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-wt-1/outs/raw_feature_bc_matrix/"

DataName <- "WT_Re2"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-wt-2/outs/raw_feature_bc_matrix/"

DataName <- "Bif3_Re1"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-bif3-1/outs/raw_feature_bc_matrix/"

DataName <- "Bif3_Re2"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped_CellRanger/Sohyun-bif3-2/outs/raw_feature_bc_matrix/"

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

gene_count_threshold <- median_gene_count - (nfold*mad_gene_count)

## 3-2) UMI cutoff
umi_count_threshold <- quantile(obj_AfterKnee$nCount_RNA, 0.10)

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




combined_plot <- grid.arrange(
  KneePlot,p1, p2, p3, p4,
  ncol = 5, 
  top = textGrob(paste0("Cell Number # after GeneCount, UMI Count, Organelle Ratio Filtering:",ncol(obj_filtered)),
                 gp = gpar(fontsize = 13))
)
# Save the combined plot
ggsave(filename = paste0("QCPlot_", DataName, ".pdf"),
       plot = combined_plot, width = 15, height = 4)


####################################
# 4) Doublet 
####################################

obj_filtered <- NormalizeData(obj_filtered)
obj_filtered <- FindVariableFeatures(obj_filtered, selection.method = "vst", nfeatures = 2000)
obj_filtered <- ScaleData(obj_filtered)
obj_filtered <- RunPCA(obj_filtered)

obj_filtered <- FindNeighbors(obj_filtered, dims = 1:10)
obj_filtered <- FindClusters(obj_filtered, resolution = 0.5)

obj_filtered <- RunUMAP(obj_filtered, dims = 1:10)

DimPlot(obj_filtered, reduction = "umap")

sweep.obj <- paramSweep(obj_filtered, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.obj, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
## knee plot
# Extract nCount_RNA values
#nCount_RNA_values <- obj_Re1$nCount_RNA

# Rank the nCount_RNA values in descending order
#ranked_nCount_RNA <- sort(nCount_RNA_values, decreasing = TRUE)

# Create a data frame with the ranks and corresponding nCount_RNA values
#knee_plot_data <- data.frame(
#  Rank = seq_along(ranked_nCount_RNA),
#  nCount_RNA = ranked_nCount_RNA
#)

# Create the knee plot
#knee_plot <- ggplot(knee_plot_data, aes(x = log10(Rank), y = log10(nCount_RNA))) +
#  geom_line() +
#  labs(
#    title = "Knee Plot",
#    x = "Rank (log10)",
#    y = "Total RNA Count (nCount_RNA) (log10)"
#  ) +
#  theme_minimal()

# Display the plot
#ggsave(filename = "kneePlot_wt_re1.pdf", plot = knee_plot, width = 4, height = 4)
