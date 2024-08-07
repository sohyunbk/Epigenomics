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

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name"),
  make_option(c("--HarmonyRDS"), type="character",
              help="HarmonyRDS", metavar="character"),
  make_option(c("--Markers"), type="character",
              help="Markers", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Harmony <- readRDS(opt$HarmonyRDS)
Name <- opt$Name

setwd(opt$WD)
Markers <- read.table(opt$Markers,head=TRUE)
#"/scratch/sb14489/3.scATAC/0.Data/MarkerGene/SeletedMarkergeneForDotPlot_RemoveUnknown.txt",head=TRUE)


generate_graph_scaleBar <- function(combined_df, GeneName, GeneID) {
  cols <- colorRampPalette(c("grey","grey","#DD3497", "#630231"), interpolate="linear", bias = .5)(100)
  
  # Define upper limit for imputed values
  upper.lim <- quantile(combined_df$value, .98, na.rm = TRUE) + 1e-6
  # Mutate the dataframe to adjust imputed values
  PlotTable <- combined_df %>%
    dplyr::mutate(final = case_when(value >= upper.lim ~ upper.lim,
                                    TRUE ~ value))
  tail(PlotTable)
  # Calculate min and max adjusted values for color scaling

  grad <- scale_colour_gradientn(colours = cols, 
                                 limits = c(0, max(PlotTable$final)))
  
  # Generate the plot
  Plot <- PlotTable  %>% 
    arrange(final)  %>% 
    ggplot(., aes(umap_1, umap_2, colour = final)) + 
    geom_jitter(position = "jitter", size = .5, stroke = 0, shape = 16) +
    scale_fill_continuous(type="viridis") + 
    grad + 
    theme_classic() + 
    ggtitle(paste0(GeneName, "\n", GeneID))

    return(Plot)
}

generate_graph <- function(combined_df, GeneName, GeneID) {
  cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
  print("Calculating Quantile")
  
  upper.lim <- quantile(combined_df$value, .98, na.rm = TRUE) + 1e-6
  
  PlotTable <- combined_df  %>%
    dplyr::mutate(final = case_when(value >= upper.lim ~ upper.lim,
                                       TRUE ~ value))
  
  min.acv <- min(PlotTable$final) - (1e-6*min(PlotTable$final))
  max.acv <- max(PlotTable$final) + (1e-6*max(PlotTable$final))
  colvec <- cols[cut(PlotTable$final, breaks=seq(min.acv, max.acv, length.out=101))]
  
  PlotTable$colors <- colvec
  
  print("Generating Graph for Marker")
  Plot <- PlotTable %>%
    arrange(final) %>%
    ggplot(aes(umap_1, umap_2, colour = colors)) +
    geom_jitter(position = "jitter", size = 1, stroke = 0, shape = 16) +
    scale_color_identity() +  # Use the actual colors in the 'colors' column
    theme_classic() +
    ggtitle(paste0(GeneName, "\n", GeneID))
  #ggsave("my_plot.png", plot = arranged_imputation_plot, width = 10, height = 8, dpi = 300)
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

#ggsave(plot=Plot,paste0(Name,"_Test.pdf"), width = 5, height = 5)
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
    gathered_list[[i]] <- generate_graph_scaleBar(combined_df, GeneName, GeneID)
    NoImputeSmooth_list[[i]] <- generate_graph_scaleBar(combined_df_Nor, GeneName, GeneID)
  }
  else{
    print("all value zero... for the gene:")
    print(GeneName)
  }
}

#ggsave(plot=Plot,paste0(Name,"_Test.pdf"), width = 5, height = 5)
captured_final_plot <- plot_grid(plotlist = gathered_list, ncol = 6)
width_cal <- 6 * 4
length_cal <- (nrow(Markers)/6 * 5)

output_name <- str_c(Name,"MarkerGeneImputed.ScaleBar", "pdf", sep = ".")
ggsave(output_name, plot = captured_final_plot,
       width = width_cal, height = length_cal,
       units = c('in'), limitsize = FALSE,
       dpi = 300)


final_plot <- plot_grid(plotlist = NoImputeSmooth_list, ncol = 6)
width_cal <- 6 * 4
length_cal <- (nrow(Markers)/6 * 5)

output_name <- str_c(Name,"MarkerGeneNotImputed.ScaleBar", "pdf", sep = ".")
ggsave(output_name, plot = final_plot,
       width = width_cal, height = length_cal,
       units = c('in'), limitsize = FALSE,
       dpi = 300)

