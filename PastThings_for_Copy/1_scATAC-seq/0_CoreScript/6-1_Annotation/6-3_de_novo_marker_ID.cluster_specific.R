library(dplyr)
library(here)
library(DESeq2)
library(tidyverse)
library(rlang)

# load arguments
args <- commandArgs(trailingOnly=T)

#args    
meta_data_file <- as.character(args[1])
gene_accessability_file <- as.character(args[2])
marker_gene_file <- as.character(args[3])
all_genes_bed <- as.character(args[4])
column_name <- as.character(args[5])
species <- as.character(args[6])
output_base <- as.character(args[7])


meta_data_file <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
gene_accessability_file <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619.txt"
marker_gene_file <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/221110_EarMarker.txt"
all_genes_bed <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed"
column_name <- "Ann_v3"
species <- "maize"
output_base <- "Ref_AnnV3"


#Load Marker Genes

message("Loading Meta Data")
marker_genes <- read.table(marker_gene_file, header =TRUE)
#Load Meta Data
head(marker_genes)

message("LOading markers Genes")
loaded_meta_data <- read.table(meta_data_file, header = TRUE) #%>% 
  #mutate(Louvain_cluster_name = str_c("LouvainC", LouvainClusters_t, sep = "_"))
head(loaded_meta_data)

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
  
  }
    
head(all_genes)
head(loaded_sparse_matric)  

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
generate_null_sample <- function(meta_data, sparse_matrix, column_name) {
  
  var <- c(column_name)
  print(var)
  print("SubSampling for Null")
  generated_subsample <- meta_data %>% 
    group_by(get(var)) %>% 
    slice_sample(prop=.50) %>% 
    ungroup() %>%
    mutate(Louvain_cluster_name = "mixed_cell_pop") %>% 
    sample_n(nrow(.))
  
  print("Splitting Null into Two groups")
  sub_sample_group_1 <- generated_subsample %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "1", sep = "."))
  
  sub_sample_group_2 <- generated_subsample %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "2", sep = "."))
  
  combined_null <- bind_rows(sub_sample_group_1, sub_sample_group_2)
  
  subsampled_sparse <- filter(sparse_matrix, barcode %in% combined_null$cellID)
  subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
  
  subsample_combined_counts <- subsample_joined %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-(!!column_name)) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  
  subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)

  return(subsample_accessability_counts_wide)
  
}
sample_cluster <- function(meta_data, sparse_matrix, cluster, column_name){
    
    var <- c(column_name)
    sprintf("sampling based off of %s and working on cluster %s", column_name, cluster)


  subsampled_meta_data <- meta_data %>% 
    dplyr::select(cellID, !!column_name) %>% 
    dplyr::filter(get(var) == cluster) %>%
    sample_n(nrow(.))
  
  group_1 <- subsampled_meta_data %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(cluster, "1", sep = "."))
  
  group_2 <- subsampled_meta_data %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(cluster, "2", sep = "."))
  
  
  combined_clustering <- bind_rows(group_1, group_2)
  
  test_set_sparse <- filter(sparse_matrix, barcode %in% combined_clustering$cellID)
  testing_join <- left_join(test_set_sparse, combined_clustering, by = c("barcode" = "cellID"))
  
  colnames(testing_join)
  
  combined_counts <- testing_join %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-!!column_name) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  accessability_counts_wide <- pivot_wider(combined_counts, names_from = "louvain_grouping_sample", 
                                           values_from = "total_accessability",   values_fill = 0) %>% 
    rowwise %>% 
    mutate(row.sum = sum(c_across(where(is.numeric)))) #%>% 
    #filter(row.sum > 5) %>% 
    #select(-row.sum)

  return(accessability_counts_wide)
}
combine_clusters_prepare <- function(sample_cluster_sparse, null_cluster_sparse){
  
  
  combined_sample_real <- left_join(sample_cluster_sparse, null_cluster_sparse, by = c("gene_name")) %>% 
    replace(is.na(.), 0)
  
  final_count_data <- as.data.frame(combined_sample_real)
  rownames(final_count_data) <- final_count_data[,1]
  final_count_data[,1] <- NULL
    
  print("DeSeq2 Design Matrix")
  print(head(final_count_data))
  
  return(final_count_data)
  
}
generate_sample_matrix <- function(combined_sparse){
  
  colnames(combined_sparse)
  generate_sample_matrix <-colnames(combined_sparse)
  one_replace <- str_replace(generate_sample_matrix, "\\.1", "")
  two_replace <- str_replace(one_replace, "\\.2", "")
  generate_sample_matrix <- two_replace 
  
  
  test <- unlist(generate_sample_matrix)
  sample_df <- as.data.frame(test)
  colnames(sample_df) <- "sample_type"
  #rownames(sample_df) <- generate_sample_matrix
  sample_df$sample_type <- factor(sample_df$sample_type)
  
  return(sample_df)
}
run_de_seq_2 <- function(counts_matrix, sample_matrix,output_base){
    #  caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, output_base)
  #counts_matrix <- final_test
  #sample_matrix <- deseq2_sample_matrix
  
    print("Inputs to Deseq2")
    print(head(counts_matrix))
    print(head(sample_matrix))

  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = sample_matrix,
                                design = ~ sample_type)
  
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  result_vector <- resultsNames(dds)
  grab_comparison_vector <- as.character(result_vector[2])
  
  
  output_name <-str_replace(grab_comparison_vector, "sample_type", "")

  res <- results(dds, name=grab_comparison_vector)
  resOrdered <- res[order(res$padj),]
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(resOrdered, rownames = "gene_name")
  
  final_output_name <- paste0(output_base,grab_comparison_vector, ".tsv")
  readr::write_tsv(as.data.frame(resOrdered.final),file=final_output_name)
  
  return(resOrdered)

}
call_cluster_specific_genes <- function(DE_seq_2_output,cluster_name,output_base){
  #  cluster_specific_genes <- call_cluster_specific_genes(caught_output, cluster_name, output_base)
   
  result_vector <- resultsNames(DE_seq_2_output)
  grab_comparison_vector <- as.character(result_vector[2])
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(DE_seq_2_output, rownames = "gene_name")
  down_regulated <- get_downregulated(resOrdered.final)
  final_output_name_2 <- paste0(output_base, ".", cluster_name, ".",grab_comparison_vector, ".cluster_specific.tsv")
  readr::write_tsv(as.data.frame(down_regulated),file=final_output_name_2)
  
  return(resOrdered.final)
  
  
}

generate_marker_list <-function(DE_seq_2_output, marker_gene_bed, total_gene_bed, cluster, output_base){
  #   generate_marker_list(cluster_specific_genes, marker_bed, all_bed, cluster_name, output_base)
  
  #DE_seq_2_output <- cluster_specific_genes
  #marker_gene_bed <- marker_bed
  #total_gene_bed <- all_bed
  #cluster <- cluster_name
  
  cluster_specific_genes.DA.2 <- get_downregulated(DE_seq_2_output) 
  cluster_specific_genes.DA <- cluster_specific_genes.DA.2
  #head(marker_gene_bed)
  #head(cluster_specific_genes.DA)
  marker_in_de_genes <- marker_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
    left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
    filter(is.na(padj) != TRUE & padj < .1 ) %>% 
    arrange(padj) %>% 
    mutate(name = str_c(name, "pval", round(padj,4), cluster, sep = "_")) %>% 
    select(chr,start,end,geneID,name,type)
  
  
  
  `%ni%` <- Negate(`%in%`)
  clust_name <- c(cluster)
  bed_in_de_genes <- total_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
    dplyr::filter(geneID %ni% as.character(marker_in_de_genes$geneID)) %>% 
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
    select(chr,start,end,geneID,name,type)
    
    print(head(marker_in_de_genes))
    print(head(bed_in_de_genes))
  #head(bed_in_de_genes)
  #head(marker_in_de_genes)
  #rbind(marker_in_de_genes, bed_in_de_genes)
  combine_markers_new_denovo <- rbind(marker_in_de_genes, bed_in_de_genes) # it has been modified

  
  output_file = paste0(output_base,".",cluster, ".markers_de_novo.visualize.bed")
  readr::write_tsv(as.data.frame(combine_markers_new_denovo),file=output_file, col_names = TRUE)


   unmodified_bed_in_de_genes <- total_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
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
    select(chr,start,end,geneID,name,type)
   
    output_file = paste0(output_base,".",cluster, ".DE_seq2.markers_de_novo.visualize.bed")
    readr::write_tsv(as.data.frame(unmodified_bed_in_de_genes),file=output_file, col_names = TRUE)

}

wrapper_function <- function(meta_data, sparse_matrix, marker_bed, all_bed, column_name, cluster_name, output_base) {
  print(cluster_name)
  function_null <- generate_null_sample(meta_data, sparse_matrix, column_name)
  sample_cluster_test <- sample_cluster(meta_data, sparse_matrix, cluster_name, column_name)
  final_test <- combine_clusters_prepare(sample_cluster_test, function_null)
  #head(final_test)
  deseq2_sample_matrix <- generate_sample_matrix(final_test)
  caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, output_base)
  cluster_specific_genes <- call_cluster_specific_genes(caught_output, cluster_name, output_base)
  #DE_seq_2_output <- caught_output
  #head(cluster_specific_genes)
  generate_marker_list(cluster_specific_genes, marker_bed, all_bed, cluster_name, output_base)
  return(cluster_specific_genes)
  
}


# Run on All --------------------------------------------------------------


WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV3"
if (!dir.exists(paste0(WD))){
  dir.create(paste0(WD))
} else {
  print("Dir already exists!")
}
setwd(paste0(WD))


print(colnames(loaded_meta_data))
head(loaded_meta_data)
final_cluster_list <- (unique(loaded_meta_data$Ann_v3))
message("Finding DA genes for the Following Clusters:")
print(final_cluster_list)

dim(loaded_sparse_matric)
head(marker_genes)
head(all_genes)
which(marker_genes$geneID=="Zm00001eb000010")
which(all_genes$geneID=="Zm00001eb000010")
all_genes[which(all_genes$geneID=="Zm00001eb000010"),]

test <- lapply(final_cluster_list, wrapper_function, meta_data = loaded_meta_data, 
               sparse_matrix=loaded_sparse_matric,  marker_bed = marker_genes, 
               column_name = "Ann_v3", 
               all_bed = all_genes, output_base = output_base)

## Done ############
############################################

## Check ###
#############################################################################
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/MtHighCutoffAllMarker_RemoveBLSelected")
Outputfile <- read.table("Ref_Denovo_Marker_MtHighCutoffAllMarker_RemoveBLSelectedsample_type_mixed_cell_pop_vs_Bundle_Sheath.tsv",header=TRUE)
#head(Outputfile)
#head(marker_genes)
Outputfile <- Outputfile[order(Outputfile$padj),]
head(Outputfile)
tail(Outputfile)
## Attach gene Symbol to Outputfile
Marker_genes <- marker_genes
colnames(Marker_genes) <- c("chr","start","end","gene_name","name","type")
#head(inner_join(Outputfile, Marker_genes, by = "gene_name"))
#head(merge(x = Outputfile, y = Marker_genes, by = "gene_name", all.x = T))

New <- merge(x = Outputfile, y = Marker_genes, by = "gene_name", all.y = T)
New <- New[order(New$padj),]
New <- New[which(New$log2FoldChange < 0),]
New[c(1:20),]
head(New)
New$name[1:20]

New <- New[order(New$log2FoldChange),]
head(New)


### All
head(Outputfile)

New <- Outputfile[order(Outputfile$padj),]
New <- New[which(New$log2FoldChange < 0),]
New[c(1:20),]
head(New)


################# Gene 
meta_orgin <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_withSub_withName.metadata.txt"
Original <- read.table(meta_orgin, header = TRUE)
head(Original)
CC<- rownames(Original)[which(Original$Ann_V1=="Companion_cell")]
meta_data_file <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBL.REF_CELLs.metadata.txt"
loaded_meta_data <- read.table(meta_data_file, header = TRUE)
head(loaded_meta_data)
loaded_meta_data$NewColor <- "Other"
loaded_meta_data[CC,]$NewColor <- "CC"

getwd()
colorr <- c("#A6CEE3", "#33A02C")
ggplot(loaded_meta_data, aes(x=umap1, y=umap2, color=factor(NewColor))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave("CCMarking_ColorRe.pdf", width=13, height=10)	




