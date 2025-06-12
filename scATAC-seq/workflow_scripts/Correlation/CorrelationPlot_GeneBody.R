library(edgeR)
library(preprocessCore)
library("edgeR")
library(tidyverse)


#GeneBodyAFile <- GeneBodyAFileBif3
#MetaFile <- MetaFileBif3
#GeneBodyAFile <- GeneBodyAFileA619
#MetaFile <- MetaFileA619
marker_correlation <- function(MetaFile, 
                               GeneBodyAFile, 
                               cluster_name = "Ann_v3"){
  ### Pull the correct meta datafrom the Socrates Object
  raw_cpm_counts_all_genes <- read.table(GeneBodyAFile,header=TRUE)
  colnames(raw_cpm_counts_all_genes) <- c("geneID","cellID","accessability")
  #raw_cpm_counts_all_genes <- sparse_matrix
  #head(sparse_matrix)
  meta_data <- read.table(MetaFile, header=TRUE)
  head(meta_data)  
  print(head(as_tibble(meta_data)))
  
  merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))  %>%
    group_by(!!sym(cluster_name), geneID)  %>%
    summarise(counts = sum(accessability, na.rm = TRUE))
  head(merged_meta_cpm_information)

  ### Alt CPM Calc
  merged_meta_cpm_information_copied <- merged_meta_cpm_information
  catch <- merged_meta_cpm_information_copied  %>%
    group_by(!!sym(cluster_name)) %>%
    group_map(~(cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>%
    unlist()
  head(catch)
  caught_values <- as_tibble(catch)
  head(caught_values)
  see <- ungroup(merged_meta_cpm_information_copied)
  idk <- bind_cols(merged_meta_cpm_information_copied,caught_values)  
  print(head(idk))
  
  merged_meta_cpm_information_copied<- bind_cols(merged_meta_cpm_information_copied,caught_values)  %>% 
    ungroup()  %>% 
    dplyr::rename(grouped_CPM = value)  %>% 
    dplyr::group_by(!!sym(cluster_name)) %>%  
    dplyr::mutate(log_cpm  = log(grouped_CPM))
  
  head(merged_meta_cpm_information_copied)
  #Apply quantile normalization
  #dim(merged_meta_cpm_information_copied)
  merged_meta.quant_norm <- merged_meta_cpm_information_copied  %>% 
    group_by(!!sym(cluster_name))  %>% 
    group_map(~(preprocessCore::normalize.quantiles(data.matrix(.x$grouped_CPM), copy = FALSE)), .keep = TRUE)  %>% 
    unlist()
  head(merged_meta.quant_norm)
  zm.quantile_normalized <- as_tibble(merged_meta.quant_norm)  %>% 
    dplyr::rename("quant_norm_cpm" = value)
  head(zm.quantile_normalized)
  
  merged_meta_cpm_information_copied <- bind_cols(merged_meta_cpm_information_copied, zm.quantile_normalized)
  table(startsWith(merged_meta_cpm_information_copied$geneID, 'Zm'))
  
  cell_type_accessability <- merged_meta_cpm_information_copied  %>% 
    dplyr::ungroup()  %>% 
    dplyr::select(!!sym(cluster_name), geneID, quant_norm_cpm)  %>% 
    pivot_wider(names_from = !!sym(cluster_name), values_from = quant_norm_cpm, values_fill = 0)
  
  head(merged_meta_cpm_information_copied)
  cell_type_accessability_logcpm <- merged_meta_cpm_information_copied  %>% 
    dplyr::ungroup()  %>% 
    dplyr::select(!!sym(cluster_name), geneID, log_cpm)  %>% 
    pivot_wider(names_from = !!sym(cluster_name), values_from = log_cpm, values_fill = 0)
  head(cell_type_accessability_logcpm)
  tail(cell_type_accessability_logcpm)
  rowsums(cell_type_accessability_logcpm)
  
  ## Other methods 
  #cell_type_accessability_logcpm
  
  
  ## Remove zero
  #cell_type_accessability_logcpm[,-1]
  #cell_type_accessability_logcpm[c(5000:5050),]
  #cell_type_accessability_logcpm_matrix <- as.matrix(cell_type_accessability_logcpm[,-1])
  #str(cell_type_accessability_logcpm_matrix)
  #head(cell_type_accessability_logcpm_matrix)
  #tail(cell_type_accessability_logcpm_matrix)
  #rowSums(cell_type_accessability_logcpm_matrix)
  #cell_type_accessability_logcpm_matrix[cell_type_accessability_logcpm_matrix==0] <- NA
  #cell_type_accessability_logcpm_matrix<-cell_type_accessability_logcpm_matrix[complete.cases(cell_type_accessability_logcpm_matrix),]
  #Cor3 <- cor(cell_type_accessability_logcpm_matrix, cell_type_accessability_logcpm_matrix,  method = "spearman")

  return(cell_type_accessability)
}

GeneBodyAFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619.txt" 
MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
cell_type_accessability_A619 <-cell_type_accessability_logcpm
  marker_correlation(GeneBodyAFileA619,MetaFileA619)

GeneBodyAFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_bif3.txt" 
MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
cell_type_accessability_bif3 <- marker_correlation(GeneBodyAFileBif3,MetaFileBif3)

cell_type_accessability_A619_matrix <-as.matrix(map(cell_type_accessability_A619[,-1], Matrix::Matrix, sparse = T)%>% 
  reduce(cbind2))
head(cell_type_accessability_A619_matrix)
colnames(cell_type_accessability_A619_matrix) <- colnames(cell_type_accessability_A619)[-1]
str(cell_type_accessability_A619_matrix)

cell_type_accessability_bif3_matrix <-as.matrix(map(cell_type_accessability_bif3[,-1], Matrix::Matrix, sparse = T)%>% 
  reduce(cbind2))
colnames(cell_type_accessability_bif3_matrix) <- colnames(cell_type_accessability_bif3)[-1]
head(cell_type_accessability_A619_matrix)
tail(cell_type_accessability_bif3_matrix)

cell_type_accessability_bif3_matrix[cell_type_accessability_bif3_matrix == "."] <- 0
cell_type_accessability_A619_matrix[cell_type_accessability_A619_matrix == "."] <- 0

rownames(cell_type_accessability_bif3_matrix) <- cell_type_accessability_bif3$geneID
rownames(cell_type_accessability_A619_matrix) <- cell_type_accessability_A619$geneID
head(cell_type_accessability_A619_matrix)
Input <- data.frame(cell_type_accessability_A619_matrix)
head(Input)
ggplot(Input, aes(x=CalloseRelated)) + 
  geom_density()
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/6.Heatmap")
ggsave("DensityPlot_CalloseRelated_logCPM.pdf", width=10, height=10)


CommonGene <- intersect(cell_type_accessability_bif3$geneID,cell_type_accessability_A619$geneID)
length(CommonGene)
dim(cell_type_accessability_A619_matrix)
dim(cell_type_accessability_bif3_matrix)
head(CommonGene)

Selected_A619 <- cell_type_accessability_A619_matrix[CommonGene,]
head(Selected_A619)
Selected_bif3 <- cell_type_accessability_bif3_matrix[CommonGene,]

Correlation <- cor(Selected_A619,Selected_bif3,  method = "pearson")
library(reshape2)
#melted_Correlation <- melt(Correlation)
#head(melted_Correlation)
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(Correlation)
melted_Correlation <- melt(upper_tri, na.rm = TRUE)
melted_Correlation <- melt(Correlation, na.rm = TRUE)

#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
ggplot(data = melted_Correlation, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "#2de361", high = "red", mid =  "#e3e32d",
                       midpoint = mean(melted_Correlation$value), 
                       limit = c(min(melted_Correlation$value),max(melted_Correlation$value)), space = "Lab", 
                       name="Pearson\nCorrelation")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() +  xlab("A619") + ylab("bif3")
#cor(x, y)
#hp       drat         wt
#mpg  -0.7761684  0.6811719 -0.8676594
#cyl   0.8324475 -0.6999381  0.7824958
#disp  0.7909486 -0.7102139  0.8879799
#> head(x)
#mpg cyl disp
#Mazda RX4         21.0   6  160
#Mazda RX4 Wag     21.0   6  160
#> head(y)
#hp drat    wt
#Mazda RX4         110 3.90 2.620
#> melted_Correlation
#Var1 Var2      value
#1  mpg   hp -0.7761684
#2  cyl   hp  0.8324475
#3 disp   hp  0.7909486
#4  mpg drat  0.6811719
#5  cyl drat -0.6999381
#6 disp drat -0.7102139
#7  mpg   wt -0.8676594
#8  cyl   wt  0.7824958
#9 disp   wt  0.8879799

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/6.Heatmap")
ggsave("A619_bif3_Ann_v3_Heatmap_notri.pdf", width=10, height=10)

library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col <- colorRampPalette(c("#2de361","#e3e32d","red"))(256)
pdf("A619_bif3_Ann_v3_Heatmap_branch.pdf", width=10, height=14)
heatmap(Correlation, scale = "none", col =  col)
dev.off()
             
