library(dplyr)
library(DESeq2)
library(tidyverse)
library(rlang)
library("optparse")

## v2 From Pablo 

# load arguments
args <- commandArgs(trailingOnly=T)

#args    
option_list = list(
  make_option(c("--meta"), type="character",
              help="meta", metavar="character"),
  make_option(c("--GeneBA"), type="character",
              help="GeneBA", metavar="character"),
  make_option(c("--marker"), type="character",
              help="marker", metavar="character"),
  make_option(c("--bed"), type="character",
              help="bed", metavar="character"),
  make_option(c("--Ann_ColumnName"), type="character",
              help="Ann_ColumnName", metavar="character"),
  make_option(c("--Species"), type="character",
              help="Species", metavar="character"),
  make_option(c("--OutputBaseName"), type="character",
              help="OutputBaseName", metavar="character"),
  make_option(c("--OutputPath"), type="character",
              help="OutputPath", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
meta_data_file <- opt$meta
gene_accessability_file <- opt$GeneBA
marker_gene_file <- opt$marker
all_genes_bed <- opt$bed
column_name_inp <- opt$Ann_ColumnName
species <- opt$Species
output_base <- opt$OutputBaseName
output_location <- opt$OutputPath

#meta_data_file_wt <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV4_metadata.txt"
#loaded_meta_data_wt <- read.table(meta_data_file_wt, header = TRUE) #%>% 

#meta_data_file <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV4_metadata.txt"
#gene_accessability_file <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_Bif3_includingZmCLE7Extended500bp.txt"
#marker_gene_file <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt"
#all_genes_bed <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr_AddZmCLE7.bed"
#column_name_inp <- "Ann_v4"
#species <- "maize"
#output_base <- "Bif3_v4"
#output_location <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV4/Bif3"

if (!dir.exists(paste0(output_location))){
  dir.create(paste0(output_location))
} else {
  print("Dir already exists!")
}
setwd(paste0(output_location))

  
#Local Development Files
#data_path <- "/Users/feilab/Projects/05.comparative_single_cell/00.data/test_data_log2fc"
#meta_data_file <- here(data_path, "sb_leaf_nmf_compressed_markers.meta.txt")
#gene_accessability_file <- here(data_path, "sorghum_bicolor.normalized_gene_acc_scores.leaf_nmf.sctGBAcounts.sparse")
#1marker_gene_file <- here(data_path, "Sb_leaf.maize_markers.ortho.visualize.bed")


#Load Marker Genes

message("Loading Meta Data")
marker_genes <- read.table(marker_gene_file, header =TRUE)
#Load Meta Data

message("LOading markers Genes")
loaded_meta_data <- read.table(meta_data_file, header = TRUE) #%>% 
#mutate(Louvain_cluster_name = str_c("LouvainC", LouvainClusters_t, sep = "_"))


if (species == "sorghum_bicolor") {
  header_names <- c("chr", "start", "end", "geneID", "val", "strand")
  #all_gene_File <- here(data_path, "Sbicolorv5.1.gene.bed")
  all_genes <- read.table(all_genes_bed, header =FALSE) %>% 
    dplyr::select(V1:V6) %>% 
    setNames(header_names) %>%  
    mutate(geneID = str_remove_all(geneID,".g")) 
  #Load Sparse Matrix
  loaded_sparse_matric <- read_delim(gene_accessability_file, delim='\t', col_names = c("gene_name", "barcode", "accessability")) %>% 
    mutate(gene_name = str_remove_all(gene_name,".g"))
  
} else {
  
  header_names <- c("chr", "start", "end", "geneID", "val", "strand")
  #all_gene_File <- here(data_path, "Sbicolorv5.1.gene.bed")
  all_genes <- read.table(all_genes_bed, header =FALSE) %>% 
    dplyr::select(V1:V6) %>% 
    setNames(header_names)
  #Load Sparse Matrix
  loaded_sparse_matric <- read_delim(gene_accessability_file, delim='\t', col_names = c("gene_name", "barcode", "accessability"))
  # https://stackoverflow.com/questions/43917904/how-to-change-column-data-type-of-a-tibble-with-least-typing
  #loaded_sparse_matric_test <- loaded_sparse_matric %<>% mutate_at(3, as.integer)
}

tail(loaded_sparse_matric)

Temp <- as.data.frame(loaded_sparse_matric,use.names=TRUE)
Temp <- Temp[-1,]
head(Temp)
str(Temp)
names(Temp)
Temp$accessability <- as.integer(Temp$accessability)
loaded_sparse_matric <- as_tibble(Temp)
head(loaded_sparse_matric)
str(loaded_sparse_matric)

# Generate Function  ------------------------------------------------------
get_downregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-.1)], rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}
get_upregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange>=.1)], rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}
generate_null_sample <- function(meta_data, sparse_matrix, slot_name, subsampled_cluster) {
  #slot_name <- "1_NA"
  slot_var <- c(slot_name)
  
  subsampled_cluster_metrics <- meta_data  %>% 
    filter(!!sym(slot_var) == subsampled_cluster)  %>% 
    ungroup()  %>% 
    dplyr::summarise(total_cells = n(), 
                     total_tn5 = sum(total))
  
  
  prcnt_10 <- subsampled_cluster_metrics$total_tn5 * .1
  up_range <- subsampled_cluster_metrics$total_tn5 + prcnt_10
  down_range <- subsampled_cluster_metrics$total_tn5 - prcnt_10
  
  passing_null_sample <- FALSE
  if (passing_null_sample == FALSE) {
    
    print("generating Null Sample")
    generated_subsample <- meta_data %>% 
      filter(!!sym(slot_var) != subsampled_cluster)  %>% 
      sample_n(nrow(.))  %>% 
      sample_n(subsampled_cluster_metrics$total_cells) %>% 
      mutate(sampled_slot = "ZZZZ_cell_pop") %>% 
      sample_n(nrow(.))
    
    generated_subsample_tn5 <- generated_subsample %>% 
      ungroup()  %>% 
      dplyr::summarise(total_tn5 = sum(total))
    
    if (between(generated_subsample_tn5, down_range, up_range) == TRUE){
      passing_null_sample <- TRUE
    } else {
      passing_null_sample <- FALSE
    }
    
  } else {
    print("Null Sample Generated")
  }
  
  
  
  
  sub_sample_group_1 <- generated_subsample %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, "QQQQQ", sep = "."))
  
  sub_sample_group_2 <- generated_subsample %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, "CCCCC", sep = "."))
  
  combined_null <- bind_rows(sub_sample_group_1, sub_sample_group_2)
  
  subsampled_sparse <- filter(sparse_matrix, barcode %in% combined_null$cellID)
  subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
  
  subsample_combined_counts <- subsample_joined %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-sampled_slot) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  
  subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)
  
  return(subsample_accessability_counts_wide)
  
}

sample_cluster <- function(meta_data, sparse_matrix, slot_name, cluster){
  
  slot_var <- c(slot_name)
  
  subsampled_meta_data <- meta_data %>% 
    dplyr::select(cellID, !!sym(slot_var)) %>% 
    dplyr::filter(!!sym(slot_var) == cluster) %>%
    sample_n(nrow(.))
  
  group_1 <- subsampled_meta_data %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), "QQQQQ", sep = "."))
  
  group_2 <- subsampled_meta_data %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), "CCCCC", sep = "."))
  
  
  combined_clustering <- bind_rows(group_1, group_2)
  
  test_set_sparse <- dplyr::filter(sparse_matrix, barcode %in% combined_clustering$cellID)
  testing_join <- left_join(test_set_sparse, combined_clustering, by = c("barcode" = "cellID"))
  
  colnames(testing_join)
  
  combined_counts <- testing_join %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-!!slot_name) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  accessability_counts_wide <- pivot_wider(combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0) %>% 
    rowwise %>% 
    mutate(row.sum = sum(c_across(where(is.numeric)))) %>% 
    filter(row.sum > 5) %>% 
    select(-row.sum)
  
  return(accessability_counts_wide)
}
combine_clusters_prepare <- function(sample_cluster_sparse, null_cluster_sparse){
  
  
  combined_sample_real <- left_join(sample_cluster_sparse, null_cluster_sparse, by = c("gene_name")) %>% 
    replace(is.na(.), 0)
  
  final_count_data <- as.data.frame(combined_sample_real)
  rownames(final_count_data) <- final_count_data[,1]
  final_count_data[,1] <- NULL
  
  return(final_count_data)
  
}
generate_sample_matrix <- function(combined_sparse){
  
  colnames(combined_sparse)
  generate_sample_matrix <-colnames(combined_sparse)
  one_replace <- str_replace(generate_sample_matrix, "\\.QQQQQ", "")
  two_replace <- str_replace(one_replace, "\\.CCCCC", "")
  generate_sample_matrix <- two_replace 
  
  
  test <- unlist(generate_sample_matrix)
  sample_df <- as.data.frame(test)
  colnames(sample_df) <- "sample_type"
  
  take_reference_factor <- "ZZZZ_cell_pop"
  alternative_factor <- as.character(unique(sample_df[sample_df$sample_type != "ZZZZ_cell_pop",]))
  
  sample_df$sample_type <- factor(sample_df$sample_type, levels = c(take_reference_factor, alternative_factor))
  
  print(sample_df)
  
  return(sample_df)
}
run_de_seq_2 <- function(counts_matrix, sample_matrix,output_base){
  
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = sample_matrix,
                                design = ~ sample_type)
  
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  result_vector <- resultsNames(dds)
  grab_comparison_vector <- as.character(result_vector[2])
  
  
  output_name <-str_replace(grab_comparison_vector, "vs_ZZZZ_cell_pop", "deseq_2_results")
  output_name <- str_replace(output_name, "sample_type_", "")
  
  res <- results(dds, name=grab_comparison_vector)
  resOrdered <- res[order(res$padj),]
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(resOrdered, rownames = "gene_name")
  
  final_output_name <- paste0(output_base, ".",output_name, ".tsv")
  readr::write_tsv(as.data.frame(resOrdered.final),file=final_output_name)
  
  return(resOrdered)
  
}
call_cluster_specific_genes <- function(DE_seq_2_output,cluster_name,output_base){
  
  result_vector <- resultsNames(DE_seq_2_output)
  grab_comparison_vector <- as.character(result_vector[2])
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(DE_seq_2_output, rownames = "gene_name")
  down_regulated <- get_upregulated(resOrdered.final)
  final_output_name_2 <- paste0(output_base, ".", cluster_name, ".upregulated_genes.deseq2_output.tsv")
  readr::write_tsv(as.data.frame(down_regulated),file=final_output_name_2)
  
  return(resOrdered.final)
}
generate_marker_list <-function(DE_seq_2_output, marker_gene_bed, total_gene_bed, cluster, output_base){
  #DE_seq_2_output <- cluster_specific_genes
  #cluster <- cluster_name
  #marker_gene_bed <- marker_bed
  #total_gene_bed <- all_bed
  cluster_specific_genes.DA.2 <- get_upregulated(DE_seq_2_output) 
  cluster_specific_genes.DA <- cluster_specific_genes.DA.2
  `%ni%` <- Negate(`%in%`)
  
  clust_name <- c(cluster)
  bed_in_de_genes <- total_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
    #dplyr::filter(geneID %ni% as.character(marker_in_de_genes$geneID)) %>% 
    left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
    filter(is.na(padj) != TRUE) %>% 
    filter(padj < .1) %>% 
    select(chr:strand,padj) %>% 
    select(-val) %>% 
    arrange(padj) %>% 
    mutate(count = row_number()) %>% 
    mutate(name = str_c(geneID, "pval", round(padj,4), clust_name, sep = "_")) %>% 
    mutate(type = clust_name) %>%
    dplyr::select(-count, -strand, -padj) %>% 
    select(chr:geneID,name,type)
  
  
  combine_markers_new_denovo <- bind_rows(bed_in_de_genes)
  
  output_file = paste0(output_base,".",cluster, ".markers_de_novo.visualize.bed")
  readr::write_tsv(as.data.frame(combine_markers_new_denovo),file=output_file, col_names = TRUE)
 
}

grab_up_regulated_genes <- function(deseq2_output, cluster_name, output_base){
  
  result_vector <- resultsNames(deseq2_output)
  grab_comparison_vector <- as.character(result_vector[2])
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(deseq2_output, rownames = "gene_name")
  cluster_specific <- get_upregulated(resOrdered.final) %>% 
    dplyr::mutate(origin = cluster_name)
  
  generate_output_name <- paste0(output_base, ".", cluster_name, ".DA_list.txt")
  
  write_delim(cluster_specific, generate_output_name, 
              col_names = TRUE, quote = "none", delim = "\t")
  
  return(cluster_specific)
}

## Replicate Aware Functions ## 

sample_cluster_replicate_aware <- function(meta_data, sparse_matrix, slot_name, cluster){
  #slot_name <- "Ann_v3"
  #cluster <- "1_NA"
  slot_var <- c(slot_name)
  subsampled_meta_data <- meta_data %>% 
    dplyr::select(cellID, library, !!sym(slot_var)) %>% 
    dplyr::filter(!!sym(slot_var) == cluster) %>%
    group_by(library) %>% 
    group_split()
  
  
  ## Shuffle the data so we can get a random sample of each
  rep_1_shuffled <- subsampled_meta_data[[1]] %>% 
    sample_n(nrow(.)) 
  
  rep_2_shuffled <- subsampled_meta_data[[2]] %>% 
    sample_n(nrow(.)) 
  
  ## SPlit replicates into two groups each (4 in all)
  rep_1_group_1 <- rep_1_shuffled %>% 
    sample_n(nrow(.)) %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), library, "QQQQQ", sep = "."))
  
  rep_1_group_2 <- rep_1_shuffled  %>% 
    sample_n(nrow(.)) %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), library, "CCCCC", sep = "."))
  
  
  rep_2_group_1 <- rep_2_shuffled %>% 
    sample_n(nrow(.)) %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), library, "QQQQQ", sep = "."))
  
  
  rep_2_group_2 <- rep_2_shuffled %>% 
    sample_n(nrow(.)) %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), library, "CCCCC", sep = "."))
  
  
  ## Bind all data together
  combined_clustering <- bind_rows(rep_1_group_1, rep_1_group_2, rep_2_group_1, rep_2_group_2)
  
  ## Join the cell IDs to the sparse matrix (actual Tn5 Counts)
  test_set_sparse <- dplyr::filter(sparse_matrix, barcode %in% combined_clustering$cellID)
  testing_join <- left_join(test_set_sparse, combined_clustering, by = c("barcode" = "cellID"))
  #tail(testing_join)
  ## Remove barcode info
  combined_counts <- testing_join %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-!!slot_name) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  ## Generate into four columns for each rep and group with counts.
  accessability_counts_wide <- pivot_wider(combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0) %>% 
    rowwise %>% 
    mutate(row.sum = sum(c_across(where(is.numeric)))) %>% 
    filter(row.sum > 2) %>% 
    select(-row.sum)
  
  return(accessability_counts_wide)
}


generate_null_sample_replicate_aware <- function(meta_data, sparse_matrix, slot_name, subsampled_cluster) {
  #slot_name <- column_name
  #head(sparse_matrix)
  #subsampled_cluster <- cluster_name
  slot_var <- c(slot_name)
  
  subsampled_cluster_metrics <- meta_data  %>% 
    filter(!!sym(slot_var) == subsampled_cluster)  %>% 
    ungroup()  %>% 
    dplyr::summarise(total_cells = n(), 
                     total_tn5 = sum(total))
  
  
  prcnt_10 <- subsampled_cluster_metrics$total_tn5 * .1
  up_range <- subsampled_cluster_metrics$total_tn5 + prcnt_10
  down_range <- subsampled_cluster_metrics$total_tn5 - prcnt_10
  
  passing_null_sample <- FALSE
  if (passing_null_sample == FALSE) {
    
    print("generating Null Sample")
    generated_subsample <- meta_data %>% 
      filter(!!sym(slot_var) != subsampled_cluster)  %>% 
      sample_n(nrow(.))  %>% 
      sample_n(subsampled_cluster_metrics$total_cells) %>% 
      mutate(sampled_slot = "ZZZZ_cell_pop") %>% 
      sample_n(nrow(.))
    
    generated_subsample_tn5 <- generated_subsample %>% 
      ungroup()  %>% 
      dplyr::summarise(total_tn5 = sum(total))
    
    result <- generated_subsample_tn5$total_tn5 >= down_range & generated_subsample_tn5$total_tn5 <= up_range
    #> down_range
    #[1] 22107964
    #> up_range
    #[1] 27020844
    #> generated_subsample_tn5
    #total_tn5
    #  20911875
    if (result == TRUE){
      passing_null_sample <- TRUE
    } else {
      passing_null_sample <- FALSE
    }
    
  } else {
    print("Null Sample Generated")
  }
  
  
  generated_subsample_rep_split <- generated_subsample  %>% 
    group_by(library) %>% 
    group_split()
  
  generated_subsample_rep_split_rep_1_shuffle <- generated_subsample_rep_split[[1]]  %>% 
    sample_n(nrow(.))
  
  generated_subsample_rep_split_rep_2_shuffle <- generated_subsample_rep_split[[2]]  %>% 
    sample_n(nrow(.))
  
  
  rep_1_sub_sample_group_1 <- generated_subsample_rep_split_rep_1_shuffle %>% 
    sample_n(nrow(.))  %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, library, "QQQQQ", sep = "."))
  
  rep_1_sub_sample_group_2 <- generated_subsample_rep_split_rep_1_shuffle %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, library, "CCCCC", sep = "."))
  
  
  rep_2_sub_sample_group_1 <- generated_subsample_rep_split_rep_2_shuffle %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, library, "QQQQQ", sep = "."))
  
  rep_2_sub_sample_group_2 <- generated_subsample_rep_split_rep_2_shuffle %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, library, "CCCCC", sep = "."))
  
  combined_null <- bind_rows(rep_1_sub_sample_group_1, rep_1_sub_sample_group_2,
                             rep_2_sub_sample_group_1, rep_2_sub_sample_group_2)
  
  subsampled_sparse <- filter(sparse_matrix, barcode %in% combined_null$cellID)
  subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
  
  subsample_combined_counts <- subsample_joined %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-sampled_slot) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)
  
  return(subsample_accessability_counts_wide)
  
}

combine_clusters_prepare_replicate_aware <- function(sample_cluster_sparse, null_cluster_sparse){
  
  
  combined_sample_real <- left_join(sample_cluster_sparse, null_cluster_sparse, by = c("gene_name")) %>% 
    replace(is.na(.), 0)
  
  final_count_data <- as.data.frame(combined_sample_real)
  rownames(final_count_data) <- final_count_data[,1]
  final_count_data[,1] <- NULL
  
  return(final_count_data)
  
}
generate_sample_matrix_replicate_aware <- function(combined_sparse){
  #combined_sparse <- final_test
  colnames(combined_sparse)
  generate_sample_matrix <-colnames(combined_sparse)
  one_replace <- str_replace(generate_sample_matrix, "\\.QQQQQ", "")
  two_replace <- str_replace(one_replace, "\\.CCCCC", "")
  generate_sample_matrix <- two_replace 
  
  
  test <- unlist(generate_sample_matrix)
  sample_df <- as.data.frame(test)
  colnames(sample_df) <- c("sample_type")
  
  
  
  final_sample_df <- as.data.frame(sample_df)  %>% 
    separate(sample_type, into = c("sample_type", "replicate"), sep ="[.]")
  
  #Have to set the factors to ensure the direction of the DE-seq2 call is correct
  take_reference_factor <- "ZZZZ_cell_pop"  
  alternative_factor_col_isoalte <- final_sample_df  %>% 
    dplyr::filter(sample_type != take_reference_factor)
  alternative_factor <- as.character(unique(alternative_factor_col_isoalte$sample_type))
  
  #Set the factors
  final_sample_df$sample_type <- factor(final_sample_df$sample_type, levels = c(take_reference_factor, alternative_factor))
  final_sample_df$replicate <- factor(final_sample_df$replicate, levels = c("bif3_Re1", "bif3_Re2"))
  
  print(final_sample_df)
  return(final_sample_df)
}

run_de_seq_2_rep_aware <- function(counts_matrix, sample_matrix, output_base){
  #run_de_seq_2_rep_aware(final_test, deseq2_sample_matrix, output_base)    
  #counts_matrix <- final_test
  #sample_matrix <- deseq2_sample_matrix
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = sample_matrix,
                                design = ~ sample_type + replicate)
  
  
  dds <- DESeq(dds)
  #head(dds)
  res <- results(dds)
  #head(res)
  #head(counts_matrix)
  #Combined <- cbind(res,counts_matrix)
  #head(Combined[which(Combined$log2FoldChange > 3 &Combined$padj < 0.05),])  
  result_vector <- resultsNames(dds)
  grab_comparison_vector <- as.character(result_vector[2])
  
  
  output_name <-str_replace(grab_comparison_vector, "vs_ZZZZ_cell_pop", "deseq_2_results")
  output_name <- str_replace(output_name, "sample_type_", "")
  
  res <- results(dds, name=grab_comparison_vector)
  resOrdered <- res[order(res$padj),]
  #head(res)  
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(resOrdered, rownames = "gene_name")
  
  final_output_name <- paste0(output_base, ".",output_name, ".tsv")
  readr::write_tsv(as.data.frame(resOrdered.final),file=final_output_name)
  
  return(resOrdered)
  
}



wrapper_function_updated <- function(meta_data, sparse_matrix, column_name,
                                     cluster_name,marker_bed, all_bed, 
                                     acr_file, output_base, output_location) {
  #cluster_name <- "CalloseRelated"
  meta_data = loaded_meta_data
  sparse_matrix=loaded_sparse_matric
  marker_bed = marker_genes
  column_name = column_name_inp
  all_bed = all_genes
  output_base = output_base
  message(paste("Working on cluster: ",cluster_name))
  #This is a gross method, but the only way to reference
  #A variable correctly when passing it as an arg
  annotation_col_quick_ref <- c(column_name)
  
  count_number_cells <- meta_data  %>% 
    group_by(!!sym(annotation_col_quick_ref))  %>% 
    summarise(cell_counts = n())  %>% 
    dplyr::filter(!!sym(annotation_col_quick_ref) == cluster_name)
  head(count_number_cells)
  
  if (count_number_cells$cell_counts < 200) {
    sample_cluster_test <- sample_cluster(meta_data, sparse_matrix, column_name, cluster_name)
    function_null <- generate_null_sample(meta_data, sparse_matrix, column_name,cluster_name)
    final_test <- combine_clusters_prepare(sample_cluster_test, function_null)
    deseq2_sample_matrix <- generate_sample_matrix(final_test)
    caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, output_base)
    
    
  } else if (count_number_cells$cell_counts > 200) {
    #Given enough cells run a replicate
    #colnames(sample_cluster_test)
    sample_cluster_test <- sample_cluster_replicate_aware(meta_data, sparse_matrix, column_name, cluster_name)
    function_null <- generate_null_sample_replicate_aware(meta_data, sparse_matrix, column_name,cluster_name)
    final_test <- combine_clusters_prepare_replicate_aware(sample_cluster_test, function_null)
    #str(final_test)
    #colnames(final_test)
    #head(final_test)
    #deseq2_sample_matrix <- final_sample_df
    deseq2_sample_matrix <- generate_sample_matrix_replicate_aware(final_test)
    caught_output <- run_de_seq_2_rep_aware(final_test, deseq2_sample_matrix, output_base)    
    
  } else {
    sample_cluster_test <- sample_cluster(meta_data, sparse_matrix, column_name, cluster_name)
    function_null <- generate_null_sample(meta_data, sparse_matrix, column_name,cluster_name)
    final_test <- combine_clusters_prepare(sample_cluster_test, function_null)
    deseq2_sample_matrix <- generate_sample_matrix(final_test)
    caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, output_base)
  }    
  
  
  
  message("Finding Upregulated genes") 
  up_reg_genes <- grab_up_regulated_genes(caught_output,cluster_name, output_base)
  #DA_marker_intersect <- intersect_ACRs_with_markers(up_reg_genes, acr_file, marker_bed)
  cluster_specific_genes <- call_cluster_specific_genes(caught_output, cluster_name, output_base)
  generate_marker_list(cluster_specific_genes, marker_bed, all_bed, cluster_name, output_base)
  return(up_reg_genes)
  
}


# Run on All --------------------------------------------------------------
print(colnames(loaded_meta_data))
final_cluster_list <- (unique(loaded_meta_data[,column_name_inp]))
message("Finding DA genes for the Following Clusters:")
print(final_cluster_list)
test <- lapply(final_cluster_list, wrapper_function_updated, 
               meta_data = loaded_meta_data, sparse_matrix=loaded_sparse_matric, 
               marker_bed = marker_genes, column_name = column_name_inp, all_bed = all_genes, 
               output_base = output_base, output_location = output_location)

#column_name_inp
#meta_data, sparse_matrix, slot_name, subsampled_cluster
