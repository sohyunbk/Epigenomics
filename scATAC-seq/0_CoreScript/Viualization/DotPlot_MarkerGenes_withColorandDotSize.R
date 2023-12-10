# Assuming you have already installed the necessary packages
library(ggplot2)

# Sample data (replace with your actual data)
set.seed(123) # For reproducibility
## In RNA-seq, the dot size is persent expression and Color is average expression.
## In ATAC-seq, the dot size is persent acc and Color is  average acc or zscore.
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
library(stringr)
library("optparse")

option_list = list(
  make_option(c("--GA"), type="character",
              help="GA", metavar="character"),
  make_option(c("--Meta"), type="character",
              help="Meta", metavar="character"),
  make_option(c("--MarkerGene"), type="character",
              help="MarkerGene", metavar="character"),
  make_option(c("--OutPathandPrefix"), type="character",
              help="OutPathandPrefix", metavar="character"),
  make_option(c("--AnnSlot"), type="character",
              help="AnnSlot", metavar="character"),
  make_option(c("--CellOrdertxt"), type="character",
              help="CellOrdertxt", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input <- opt$GA
meta <- opt$Meta
gene <- opt$MarkerGene
OutputPathandName<- opt$OutPathandPrefix
slot_var <- opt$AnnSlot
CellOrder <- opt$CellOrdertxt

#input <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re.txt"
#meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt"
#gene <- "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt"
#OutputPathandName <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/7.DotPlot/A619_Annv4_DenovoGenes"
#slot_var <- "Ann_v4"
#CellOrder <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt"


message("Loading Data....")
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)
gene_markers <- gene_markers  %>%
  arrange(type)

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
                            by = c(slot_var, "geneID"))

# Print the first few rows of the combined table
head(combined_table)
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


filtered_table <- result_table %>%
  filter(grepl(paste(gene_markers$name, collapse = "|"), name, ignore.case = TRUE))
#filtered_table[,c(6:13)]

CellTypeOrder <- rev(readLines(CellOrder))
filtered_table$Ann <- factor(filtered_table[[slot_var]], levels = CellTypeOrder)

# Function to find the first matching element in CellTypeOrder
first_match <- function(name, order_list) {
  for (order_item in order_list) {
    if (str_detect(name, order_item)) {
      return(order_item)
    }
  }
  return(NA)  # Return NA if no match is found
}
# Create a new column for sorting 
## type is cell type matched with the marker gene!!
## Ann is Annotated cell type
filtered_table$sorting_key <- sapply(filtered_table$type, first_match, order_list = CellTypeOrder)
# Convert the sorting_key column to a factor with levels in the order of CellTypeOrder
filtered_table$sorting_key <- factor(filtered_table$sorting_key, levels = CellTypeOrder)
# Arrange the table by the sorting_key column
ordered_table <- filtered_table %>% 
  arrange(sorting_key, type)

ordered_table$name <- factor(ordered_table$name,levels=rev(unique(ordered_table$name)))
# Print the resulting table
ggplot(ordered_table, aes(x = name, y = Ann, 
                        size = proportion_expressing, color = Zscore[,1])) +
  geom_point() +
  #scale_size_continuous(range = c(1, 10)) +  # Adjust the range for the size of dots
  labs( x = "Gene",
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


ggsave(paste0(OutputPathandName,".pdf"), width=0.3*length(unique(filtered_table$geneID)),
       height=4)

