## Plot subset of marker genes ##

# load arguments
args <- commandArgs(T)
#if(length(args)!=5){stop("Rscript normGBA.R <gene.sparse> <meta> <Zea_mays.AGPv4.36.Allgene.nuclear.bed> <prefix> <F>")}
input <- as.character(args[1])
meta <- as.character(args[2])
gene <- as.character(args[3])
prefix <- as.character(args[4])

# load libraries
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


input <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_IncludingZmCLE7/opt_allgenes_impute.activity.rds"
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV4_metadata.txt"

input <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_IncludingZmCLE7/opt_allgenes_impute.activity.rds"
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV4_metadata.txt"

gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt"
prefix <- "A619_includingZmCLE7"

#sparse_matrix <-  obj 
sparse_matrix <- readRDS(input)
head(sparse_matrix)
imputed_sparse <- sparse_matrix
imputed_sparse <- (sparse_matrix$impute.activity)
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)  

head(imputed_sparse[,c(1:10)])
dim(imputed_sparse)
gene_markers <- gene_markers  %>% 
  arrange(type)

all_markers <- gene_markers$geneID
all_markers <- all_markers[all_markers %in% rownames(imputed_sparse)]
all_markers <- unique(all_markers[all_markers %in% rownames(imputed_sparse)])

gathered_list <- list()

take_length_markers <- seq(1,length(all_markers))
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_IncludingZmCLE7")

for ( z in take_length_markers ) {
  print(all_markers[[z]])
  print(z)
  i <- all_markers[[z]]
  
  storage_character <- as.character(z)
  
  print("Grabbing Sparse Meta Info")
  grab_gene_info <- gene_markers[gene_markers$geneID == i, ]
  pull_sparse_matrix_gene <- imputed_sparse[i,]
  
  print("Converting to Tibble")
  imputer_sparse_2_filtered <- as_tibble(pull_sparse_matrix_gene, rownames = "cellID")
  
  
  imputated_sparse_2_join <- left_join(meta_data, imputer_sparse_2_filtered, by = c("cellID"), copy = TRUE)  %>% 
    replace(is.na(.), 0)
  
  #cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
  cols <- colorRampPalette(c("grey","grey","grey","grey","grey","#DD3497", "#630231"), interpolate="linear", bias = .5)(100)
  
  print("Calculating Quantile")
  min.acv <- min(acv) - (1e-6*min(acv))
  max.acv <- max(acv) + (1e-6*max(acv))
  
  upper.lim <- quantile(imputated_sparse_2_join$value, .90, na.rm = TRUE) + 1e-6
  
  graphable_sparse_2 <- imputated_sparse_2_join  %>% 
    dplyr::mutate(final_ac = case_when(value >= upper.lim ~ upper.lim,
                                       TRUE ~ value))
  
  grad <- scale_colour_gradientn(colours = cols, 
                                 limits = c(0, max(graphable_sparse_2$final_ac)))
  
  generate_title <- str_c(grab_gene_info$name, grab_gene_info$geneID, sep ="\n")
  
  print("Generating Graph for Marker")
  arranged_imputation_plot <- graphable_sparse_2  %>% 
    arrange(final_ac)  %>% 
    ggplot(., aes(umap1, umap2, colour = final_ac)) + 
    geom_jitter(position = "jitter", size = .5, stroke = 0, shape = 16) +
    scale_fill_continuous(type="viridis") + 
    grad + 
    theme_classic() + 
    #theme(legend.position="none") +
    ggtitle(generate_title)
  
  gathered_list[[storage_character]] <- arranged_imputation_plot
  
}


print("Plotting All Markers in Grid...")
captured_final_plot <- plot_grid(plotlist = gathered_list, ncol = 6)
width_cal <- 6 * 5
length_cal <- (length(all_markers)/6 * 5)


output_name <- str_c(prefix, "pdf", sep = ".")

ggsave(output_name, plot = captured_final_plot, 
       width = width_cal, height = length_cal, 
       units = c('in'), limitsize = FALSE,
       dpi = 300)
