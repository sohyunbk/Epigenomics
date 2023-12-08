# Assuming you have already installed the necessary packages
library(ggplot2)

# Sample data (replace with your actual data)
set.seed(123) # For reproducibility

## In RNA-seq, the dot size is persent expression and Color is average expression.
## In ATAC-seq, the dot size is persent acc and Color is  average acc or zscore.

## From Pablo

# load libraries
library(dplyr)
library("edgeR")
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(tidyverse)


#working_dir <- "/scratch/jpm73279/comparative_single_cell/08.annotation_figures/zea_mays"

# load arguments
args <- commandArgs(T)
#if(length(args)!=5){stop("Rscript normGBA.R <gene.sparse> <meta> <Zea_mays.AGPv4.36.Allgene.nuclear.bed> <prefix> <F>")}
input <- as.character(args[1])
meta <- as.character(args[2])
gene <- as.character(args[3])
annot_col <- as.character(args[4])
prefix <- as.character(args[5])


slot_var <- c(annot_col)

#gene_bodysorghum_bicolor_tis_leaf_nmf_step_2_knn_27.out.de_novo.rds
input <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re.txt"
meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/221130_EarMarker.txt"
slot_var <- "Ann_v3"
#gene_DA <- here(working_dir,"00.data/Zm-B73-REFERENCE-NAM_Zm00001eb.1.genes.bed")
#prefix <- "TEST_SORGHUM_TEST"


message("Loading Data....")
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)
gene_markers <- gene_markers  %>%
  arrange(type)

all_markers <- gene_markers$geneID

#raw_cpm_counts_all_genes <- read_delim(input, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
#  dplyr::mutate(cellID = barcode)  %>%
#  dplyr::mutate(geneID = gene_name)

raw_cpm_counts_all_genes <- read_delim(input, delim="\t") %>%
  dplyr::mutate(cellID = barcode)  %>%
  dplyr::mutate(geneID = gene_name)

head(raw_cpm_counts_all_genes)
colnames(meta_data)


message("Calculating CPM Values....")
merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))  %>%
  group_by(!!sym(slot_var), geneID)  %>%
  summarise(counts = sum(accessability, na.rm = TRUE))

### Alt CPM Calc
merged_meta_cpm_information_copied <- merged_meta_cpm_information
catch <- merged_meta_cpm_information_copied  %>%
  group_by(!!sym(slot_var)) %>%
  group_map(~(cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>%
  unlist()

caught_values <- as_tibble(catch)
see <- ungroup(merged_meta_cpm_information_copied)
merged_meta_cpm_information_copied<- bind_cols(merged_meta_cpm_information_copied,caught_values)  %>% 
  rename(grouped_CPM = value)


message("Calculating ZScore Approximations...")
head(merged_meta_cpm_information_copied)
altered_deseq2 <- merged_meta_cpm_information_copied %>% 
  dplyr::select(-counts) %>% 
  pivot_wider(names_from = geneID, values_from = grouped_CPM, values_fill = 0) %>% 
  pivot_longer(cols = -!!sym(slot_var), names_to = "geneID", values_to = "grouped_CPM") %>% 
  group_by(geneID) %>% 
  mutate(Zscore = scale(grouped_CPM)) %>% 
  ungroup()  %>% 
  #mutate(relative_accessability = rescale(Zscore, to = c(0,1))) %>% 
  group_by(!!sym(slot_var))  %>% 
  mutate(Zscore_group = scale(Zscore))

message("Calculating Proportion of cells marker is Accessible IN...")
# Create Proportion Cells Accessible Metrics ------------------------------
merged_meta_cellID_values <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))
take_unq_genes <- unique(merged_meta_cellID_values$geneID)

head(merged_meta_cellID_values)

merged_meta_cellID_values_all_genes <- merged_meta_cellID_values %>% 
  select(cellID, !!sym(slot_var), accessability, geneID) 


wider_all_genes_altered <- merged_meta_cellID_values_all_genes %>% 
  pivot_wider(names_from = geneID, 
              values_from = accessability,  
              values_fill = 0) %>% 
  pivot_longer(cols = c(-!!sym(slot_var), -cellID), 
               names_to = "geneID", 
               values_to = "accessability") %>% 
  mutate(expression_bool = case_when(accessability < 1 ~ 0,
                                     accessability >= 1 ~ 1)) %>% 
  group_by(!!sym(slot_var), geneID) %>% 
  summarise(total_cells = n(), 
            proportion_expressing = (sum(expression_bool)/total_cells * 100))
wider_all_genes_altered
altered_deseq2

# Remove the grouping from the tables before combining
wider_all_genes_altered <- wider_all_genes_altered %>% ungroup()
altered_deseq2 <- altered_deseq2 %>% ungroup()

# Combine the tables without duplicating the common columns "Ann_v3" and "geneID"
combined_table <- left_join(wider_all_genes_altered, altered_deseq2, 
                            by = c("Ann_v3", "geneID"))

# Print the first few rows of the combined table
print(combined_table)
# Print the first few rows of the combined table

##################################################
### Make input:  gene - cell_type - Acc percent - zscore.
###################################################

##########
head(gene_markers)
SelectGenes <- gene_markers
# Perform a left join to add the "name" column to the combined_table based on geneID
result_table <- inner_join(combined_table, gene_markers, by = "geneID")
# Print the first few rows of the resulting table
## Narrow Down more genes!!
print(result_table,width=Inf)
genesIntrested <- c("fdl", "ocl4",  "lg2", "zag1",
                    "wus1","RA1", "ra2","OFB1", "gif1", "wox4", "ZmMP_1",
                   "vnd6", "xcp1", "ZmSMXL3", "ZmSMXL4", "nac78", "ZmNEN1",
                   "cyc3", "cyc4")


filtered_table <- result_table %>%
  filter(grepl(paste(genesIntrested, collapse = "|"), name, ignore.case = TRUE))

filtered_table$name <- factor(filtered_table$name, levels = genesIntrested)
levels(factor(filtered_table$Ann_v3))
CellTypeOrder <- rev(c("L1","L1atFloralMeristem","FloralMeristem_SuppressedBract",
                   "IM-OC","SPM-base_SM-base","IM_SPM_SM",
                   "ProcambialMeristem_ProtoXylem_MetaXylem",
                   "PhloemPrecursor", "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma",
                   "XylemParenchyma_PithParenchyma",
                   "BundleSheath_VascularSchrenchyma",
                   "CalloseRelated","G2_M"))
filtered_table$Ann_v3 <- factor(filtered_table$Ann_v3, levels = CellTypeOrder)

# Print the resulting table
print(filtered_table)

ggplot(filtered_table, aes(x = name, y = Ann_v3, 
                        size = proportion_expressing, color = Zscore[,1])) +
  geom_point() +
  #scale_size_continuous(range = c(1, 10)) +  # Adjust the range for the size of dots
  labs(title = "Dot Plot with Gene Expression and Other Data",
       x = "Gene",
       y = "Cell Type",
       size = "Percent ACC",
       color = "Zscore") +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove the grid
        axis.line = element_line(color = "black"),  # Add x-axis and y-axis lines
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # Remove the panel border
  
  # Change the color scale to range from blue to white to red
  scale_color_gradient2(low = "blue",
                        mid = "white", high = "red", 
                        midpoint = mean(filtered_table$Zscore[,1]))

ggsave("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/7.DotPlot/A619_dotplot.pdf"
       , width=10, height=4.5)

