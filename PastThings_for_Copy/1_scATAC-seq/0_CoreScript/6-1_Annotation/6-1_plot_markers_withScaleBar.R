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

#gene_bodysorghum_bicolor_tis_leaf_nmf_step_2_knn_27.out.de_novo.rds
#input <- here("/home/jpm73279/r_script_dev/lw_plotting","gene_bodysorghum_bicolor_tis_leaf_nmf_step_2_knn_26.out.rds")
#input_2 <- here("/home/jpm73279/r_script_dev/lw_plotting","gene_bodysorghum_bicolor_tis_leaf_nmf_step_3_knn_27.out.rds")
#meta <- here("/home/jpm73279/r_script_dev/lw_plotting","Sb_leaf.merged_replicates.NMF.full.metadata.txt")
#gene <- here("/home/jpm73279/r_script_dev/lw_plotting","Sb_leaf.maize_markers.ortho.visualize.bed")
#gene_DA <- here("/home/jpm73279/r_script_dev/lw_plotting","sorghum_bicolor.tis_leaf_nmf.cluster.DA_genes.merged.bed")
#prefix <- "TEST_SORGHUM_TEST"
setwd("/scratch/sb14489/3.scATAC_flo/Test_Marker_withBar")

input <- "/scratch/sb14489/3.scATAC_flo/7.CellAnnotation/Ref_500bpBin/out.rds"
meta <- "/scratch/sb14489/3.scATAC_flo/5.Socrates/4_CombineAll_AfterD/Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.REF_CELLs.ResKN_0.7_30_metadata.txt"
gene <- "/scratch/sb14489/3.scATAC_flo/7.CellAnnotation/0.MarkerGene/Ear.markers.final.txt"
prefix <- "Test_Marker_withBar"

#sparse_matrix <-  obj 
sparse_matrix <- readRDS(input)
imputed_sparse <- (sparse_matrix$impute.activity)
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)  #%>% 
    #mutate(geneID = case_when((grepl(".m", geneID) ~ str_c(geneID,"g", sep=".")), 
    #                          TRUE ~ geneID))

head(imputed_sparse[c(1:10),])
dim(imputed_sparse)
gene_markers <- gene_markers  %>% 
    arrange(type)


all_markers <- gene_markers$geneID
all_markers <- all_markers[all_markers %in% rownames(imputed_sparse)]
all_markers <- unique(all_markers[all_markers %in% rownames(imputed_sparse)])
all_markers <- all_markers[all_markers != "ta.jpa7.902.35711FPf1"]
all_markers <- all_markers[all_markers != "ta.jpa7.902.14508FPf3"]

gathered_list <- list()



take_length_markers <- seq(1,length(all_markers))


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
    cols <- colorRampPalette(c("grey80","#DD3497", "#49006A"), interpolate="linear", bias = .5)(100)
    
    print("Calculating Quantile")
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
