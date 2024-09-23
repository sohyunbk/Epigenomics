library("here")
library(devtools)
library(Seurat)
library(RANN)
load_all('/home/sb14489/Socrates')
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(sctransform)
library(tidyverse)
library(here)
library(patchwork)
library(cowplot)
library("optparse")
library(rlang)
library(ggplot2)

Harmony <- readRDS("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/4.Harmony/obj_afterHarmony_WTRe1andRe2_UMI1000.rds")
Name <- "WTRe1andRe2_UMI1000"

setwd("/scratch/sb14489/4.scRNAseq/2.snRNA-seq/5.MarkerGene")
Markers <- read.table("/scratch/sb14489/3.scATAC/0.Data/MarkerGene/SeletedMarkergeneForDotPlot_RemoveUnknown.txt",head=TRUE)
head(Markers)
dim(Markers)

MarkerGenePlot <- 
  FeaturePlot(Harmony, features = Markers$geneID)

ggsave(plot=MarkerGenePlot,paste0(Name,"_28MarkerGene.pdf"), width = 10, height = 22)
### function for figure ##

generate_graph <- function(combined_df, GeneName, GeneID) {
  # Generate color palette
  #WhichValue <- "ImputedValue"
  cols <- colorRampPalette(c("grey","grey","grey","grey","grey","#DD3497", "#630231"), interpolate="linear", bias = .5)(100)
  
  # Define upper limit for imputed values
  upper.lim <- quantile(combined_df$value, .98, na.rm = TRUE) + 1e-6
  
  # Mutate the dataframe to adjust imputed values
  graphable_sparse_2 <- combined_df %>%
    dplyr::mutate(final = case_when(value >= upper.lim ~ upper.lim,
                                    TRUE ~ value))
  
  # Calculate min and max adjusted values for color scaling
  min.acv <- min(graphable_sparse_2$final) - (1e-6 * min(graphable_sparse_2$final))
  max.acv <- max(graphable_sparse_2$final) + (1e-6 * max(graphable_sparse_2$final))
  
  # Assign colors based on the adjusted values
  colvec <- cols[cut(graphable_sparse_2$final, breaks = seq(min.acv, max.acv, length.out = 101))]
  
  # Add colors to the dataframe
  graphable_sparse_2$colors <- colvec
  
  print("Generating Graph for Marker")
  
  # Generate the plot
  Plot <- graphable_sparse_2  %>% 
    arrange(final)  %>% 
    ggplot(., aes(umap_1, umap_2, colour = final)) + 
    geom_jitter(position = "jitter", size = .5, stroke = 0, shape = 16) +
    scale_fill_continuous(type="viridis") + 
    grad + 
    theme_classic() + 
    #theme(legend.position="none") +
    ggtitle(paste0(GeneName, "\n", GeneID))
  # Return the list containing the plot
  return(Plot)
}

############################
## smooth and imputation 
############################


k=25
step=2 # originally 3
npcs=19
pca_embeddings <- Embeddings(Harmony, reduction = "pca")
# do l2-norm
pcs <- apply(pca_embeddings, 2, function(x){x/sqrt(sum(x^2))})
head(pcs)

# get KNN
message("   * finding knn 
        graph ...")
knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
j <- as.numeric(x = t(x = knn.graph))
i <- ((1:length(x = j)) - 1) %/% k + 1
edgeList = data.frame(i, j, 1)
A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])


# Smooth graph
##some questions
##1. what's meanning of A
##2. A/Matrix::rowSums(A)
##3. meanning of A%*%A
message("   * smoothing graph ...")
A = A + t(A) ##A now is symmetric
A = A / Matrix::rowSums(A)
step.size = step
if(step.size > 1){
  for(i in 1:step.size){
    message("     ~ step ",i)
    A = A %*% A
  }
}
head(A)
dim(A)
head(pca_embeddings[,c(1:10)])
colnames(A) <- rownames(pca_embeddings)
saveRDS(A,paste0(Name,".MarkovMatrix.rds"))

############################################

message("   * smoothing activity ...")
data <- GetAssayData(object = Harmony, layer = "data") # Normalized data matrix
#data <- GetAssayData(object = Harmony, layer = "counts") # raw data matrix
head(data[,c(1:10)])
out <- t(A %*% t(data))
colnames(out) <- colnames(data)
rownames(out) <- rownames(data)
print(head(out[,1:10]))
dim(out)

head(out["Zm00001eb000010",])
head(Harmony@meta.data)
umap_embeddings <- Embeddings(Harmony, reduction = "umap")
head(umap_embeddings)
str(combined_df)


#### Draw all Marker genes
gathered_list <- list()
NoImputeSmooth_list <- list()

for (i in c(1:nrow(Markers))){
  GeneID <- Markers[i,"geneID"]
  GeneName <- Markers[i,"name"]
  combined_df <- data.frame(cbind(umap_embeddings, value =out[GeneID,]))
  combined_df_Nor <- data.frame(cbind(umap_embeddings, value = data[GeneID,] ))
  head(combined_df)
  head(combined_df_Nor)
  ########################################
  if(sum(combined_df_Nor$value)!=0){
    gathered_list[[i]] <- generate_graph(combined_df, GeneName, GeneID)
    NoImputeSmooth_list[[i]] <- generate_graph(combined_df_Nor, GeneName, GeneID)
  }
  else{
    print("all value zero... for the gene:")
    print(GeneName)
  }
}

#ggsave(plot=arranged_imputation_plot,paste0(Name,"_Test.pdf"), width = 10, height = 22)
captured_final_plot <- plot_grid(plotlist = gathered_list, ncol = 6)
width_cal <- 6 * 4
length_cal <- (nrow(Markers)/6 * 5)

output_name <- str_c(Name,"MarkerGeneImputed", "pdf", sep = ".")
ggsave(output_name, plot = captured_final_plot,
       width = width_cal, height = length_cal,
       units = c('in'), limitsize = FALSE,
       dpi = 300)


final_plot <- plot_grid(plotlist = NoImputeSmooth_list, ncol = 6)
width_cal <- 6 * 4
length_cal <- (nrow(Markers)/6 * 5)

output_name <- str_c(Name,"MarkerGeneNotImputed", "pdf", sep = ".")
ggsave(output_name, plot = final_plot,
       width = width_cal, height = length_cal,
       units = c('in'), limitsize = FALSE,
       dpi = 300)

