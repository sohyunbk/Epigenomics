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
library("optparse")
library(rlang)
library(ggplot2)

#e.g.
#meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Cluster1_Recluster_Sub_res1_knear100_Partmetadata.txt"
#gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/221130_EarMarker.txt"
#OutputPath <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3SubClsters"

option_list = list(
  make_option(c("--imputed_sparse"), type="character", 
              help="imputed_sparse", metavar="character"),
  make_option(c("--meta"), type="character", 
              help="meta"),
  make_option(c("--gene"), type="character",
              help="gene", metavar="character"),
  make_option(c("--prefix"), type="character", 
              help="prefix", metavar="character"),
  make_option(c("--OutputPath"), type="character", 
              help="OutputPath", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

imputed_sparse_rds <- opt$imputed_sparse
meta <- opt$meta
gene <- opt$gene
prefix <- opt$prefix
OutputPath <- opt$OutputPath

imputed_sparse <- readRDS(imputed_sparse_rds)
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene) 

head(gene_markers)
all_markers <- gene_markers$geneID
all_markers <- all_markers[all_markers %in% rownames(imputed_sparse)]
all_markers <- unique(all_markers[all_markers %in% rownames(imputed_sparse)])
length(all_markers)

gathered_list <- list()

print(all_markers)
take_length_markers <- seq(1,length(all_markers))
#z<-1

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

    cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
    #cols <- colorRampPalette(c("grey80","#DD3497", "#49006A"), interpolate="linear", bias = .5)(100)

    print("Calculating Quantile")
    global_acc <- imputed_sparse[rownames(imputed_sparse) == i]
    upper.lim <- quantile(global_acc, .98, na.rm = TRUE) + 1e-6

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
            geom_jitter(position = "jitter", size = 1, stroke = 0, shape = 16) +
            scale_fill_continuous(type="viridis") +
            grad +
            theme_classic() +
            #theme(legend.position="none") +
            ggtitle(generate_title)

    gathered_list[[storage_character]] <- arranged_imputation_plot

}


print("Plotting All Markers in Grid...")
captured_final_plot <- plot_grid(plotlist = gathered_list, ncol = 6)
width_cal <- 6 * 4
length_cal <- (length(all_markers)/6 * 4)

setwd(OutputPath)

output_name <- str_c(prefix, "pdf", sep = ".")
ggsave(output_name, plot = captured_final_plot,
       width = width_cal, height = length_cal,
       units = c('in'), limitsize = FALSE,
       dpi = 300)


