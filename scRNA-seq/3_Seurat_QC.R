library(stringi)
library(SeuratData)
library(Seurat)
#library(Signac)
library(ggplot2)
library(cowplot)
library(tidyr)
library(gridExtra)
library(grid)  # Load the grid package

setwd("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/3.Seurat")

DataName <- "WT_Re1"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped/Sohyun-wt-1/outs/filtered_feature_bc_matrix/"

DataName <- "WT_Re2"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped/Sohyun-wt-2/outs/filtered_feature_bc_matrix/"

DataName <- "Bif3_Re1"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped/Sohyun-bif3-1/outs/filtered_feature_bc_matrix/"

DataName <- "Bif3_Re2"
InputFile <- "/scratch/sb14489/4.scRNAseq/2.snRNA-seq/2.Mapped/Sohyun-bif3-2/outs/filtered_feature_bc_matrix/"

data <- Read10X(data.dir = InputFile)
obj <- CreateSeuratObject(counts = data, project = DataName, min.cells = 100, min.features = 200)
obj[["percent.organelle"]] <- PercentageFeatureSet(obj, pattern = "^gene:GRMZM")

obj$log10_nFeature_RNA <- log10(obj$nFeature_RNA + 0.00001)  # Add 1 to avoid log(0)
obj$log10_nCount_RNA <- log10(obj$nCount_RNA + 0.00001)  # Add 1 to avoid log(0)


### 1) Gene count - feature count filter
median_gene_count <- median(obj$nFeature_RNA)
mad_gene_count <- mad(obj$nFeature_RNA)

gene_count_threshold <- median_gene_count

## 2) UMI cutoff
umi_count_threshold <- quantile(obj_filtered_genecount$nCount_RNA, 0.20)

## 3) Org cutoff 
nOrg <- 5

obj_filtered <- subset(obj, subset = nCount_RNA > umi_count_threshold &
                         nFeature_RNA > gene_count_threshold &
                         percent.organelle < nOrg
                         )




### violet Plot
p1 <- ggplot(obj@meta.data, aes(x=orig.ident, y=log10_nFeature_RNA)) + 
  geom_violin()+
  geom_hline(yintercept = log10(gene_count_threshold), linetype = "dashed", color = "red") + 
  theme_minimal() + 
  labs(x = DataName, y = "Log10 Gene Count")

p2 <- ggplot(obj@meta.data, aes(x=orig.ident, y=log10_nCount_RNA)) + 
  geom_violin()+
  geom_hline(yintercept = log10(umi_count_threshold), linetype = "dashed", color = "red") + 
  theme_minimal() + 
  labs(x = DataName, y = "Log10 Read Count")

p3 <- ggplot(obj@meta.data, aes(x=orig.ident, y=percent.organelle)) + 
  geom_violin()+
  geom_hline(yintercept = nOrg, linetype = "dashed", color = "red") + 
  theme_minimal() + 
  labs(x = DataName, y = "Organelle Percentage")

combined_plot <- grid.arrange(
  p1, p2, p3, 
  ncol = 3, 
  top = textGrob(paste0("Cell # before Filtering:",ncol(obj),"\nCell # after Filtering:",ncol(obj_filtered)),
                 gp = gpar(fontsize = 16, fontface = "bold"))
)
# Save the combined plot
ggsave(filename = paste0("violin_plot_", DataName, ".pdf"), plot = combined_plot, width = 15, height = 4)







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
