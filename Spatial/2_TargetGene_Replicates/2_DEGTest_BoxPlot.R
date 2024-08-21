library(edgeR)
library(ggplot2)


##### 1) load Data and write the Replicate info!
file_path <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMAM/"

# WT Re1
data <- read.csv(paste0(file_path,"A619_B_Cluser5.csv"), row.names = 1, header = TRUE)
#data <- as.data.table(data)
WT_Re1 <- colSums(data)

#WT Re2
data <- read.csv(paste0(file_path,"A619_C_Cluser1.csv"), row.names = 1, header = TRUE)
WT_Re2_1 <- colSums(data)
data <- read.csv(paste0(file_path,"A619_C_Cluser7.csv"), row.names = 1, header = TRUE)
WT_Re2_2 <- colSums(data)
data <- read.csv(paste0(file_path,"A619_D_WTCZOnly_new.csv"), row.names = 1, header = TRUE)
WT_Re2_3 <- colSums(data)
WT_Re2 <- WT_Re2_1+WT_Re2_2+WT_Re2_3

# Bif3 Re1
data <- read.csv(paste0(file_path,"A619_D_Bif3CZOnly.csv"), row.names = 1, header = TRUE)
Bif3_Re1 <- colSums(data)

# Bif3 Re2
data <- read.csv(paste0(file_path,"Bif3_A_Cluser1.csv"), row.names = 1, header = TRUE)
Bif3_Re2_1 <- colSums(data)
data <- read.csv(paste0(file_path,"Bif3_A_Cluser10.csv"), row.names = 1, header = TRUE)
Bif3_Re2_2 <- colSums(data)
Bif3_Re2 <- Bif3_Re2_1+Bif3_Re2_2

# Bif3 Re3
data <- read.csv(paste0(file_path,"Bif3_B_Cluser5.csv"), row.names = 1, header = TRUE)
Bif3_Re3_1 <- colSums(data)
data <- read.csv(paste0(file_path,"Bif3_B_Cluser6.csv"), row.names = 1, header = TRUE)
Bif3_Re3_2 <- colSums(data)
Bif3_Re3 <- Bif3_Re3_1+Bif3_Re3_2

# Bif3 Re4
data <- read.csv(paste0(file_path,"Bif3_C_Cluser11.csv"), row.names = 1, header = TRUE)
Bif3_Re4_1 <- colSums(data)
data <- read.csv(paste0(file_path,"Bif3_C_Cluser5.csv"), row.names = 1, header = TRUE)
Bif3_Re4_2 <- colSums(data)
data <- read.csv(paste0(file_path,"Bif3_D_Cluser10.csv"), row.names = 1, header = TRUE)
Bif3_Re4_3 <- colSums(data)
data <- read.csv(paste0(file_path,"Bif3_D_Cluser7.csv"), row.names = 1, header = TRUE)
Bif3_Re4_4 <- colSums(data)
Bif3_Re4 <- Bif3_Re4_1+Bif3_Re4_2+Bif3_Re4_3+Bif3_Re4_4

CountTable <- data.frame(WT_Re1,WT_Re2,Bif3_Re1,Bif3_Re2,Bif3_Re3,Bif3_Re4)
CountTable <- CountTable[grep("^Zm00", rownames(CountTable)), ]
head(CountTable)
tail(CountTable)

### EdgeR count info: row - gene column Sample
group <- factor(c(rep("WT",2), rep("Bif3",4)))
dge <- DGEList(counts=CountTable, group=group)
dge <- calcNormFactors(dge, method="TMM")
tmm_normalized_counts <- cpm(dge, normalized.lib.sizes=TRUE)
head(tmm_normalized_counts)

AllPlotgenes_symbol <- c('arftf4', 'arftf30', 'arftf18', 'arftf3', 'arftf20', 'arftf25', 'arftf10', 'arftf36', 'arftf23', 'arftf26', 'knox1')
AllPlotgenes <- c('Zm00001eb433460', 'Zm00001eb066640', 'Zm00001eb067270', 'Zm00001eb142540', 'Zm00001eb224680', 'Zm00001eb232120', 'Zm00001eb243930', 'Zm00001eb292830', 'Zm00001eb363810', 'Zm00001eb370810', 'Zm00001eb001720')

tmm_normalized_selected <- tmm_normalized_counts[AllPlotgenes,]

Plotlist <- list()
i <-11
i <- 9
for (i in c(1:nrow(tmm_normalized_selected))) {
  # Extract gene data for WT and Bif3
  gene <- rownames(tmm_normalized_selected)[[i]]
  gene_data <- tmm_normalized_selected[gene, ]
  wt_values <- gene_data[c(1:2)]
  bif3_values <- gene_data[c(3:6)]
  
  # Create a data frame for ggplot
  plot_data <- data.frame(
    value = c(wt_values, bif3_values),
    group = rep(c("WT", "Bif3"), c(length(wt_values), length(bif3_values)))
  )
  
  # Plotting the violin plot
  p <- ggplot(plot_data, aes(x = group, y = value)) +
    geom_boxplot() +
    labs(title = AllPlotgenes_symbol[[i]], x = "Group", y = "Expression Level") +
    theme_minimal()
  Plotlist[[i]] <- p
}

library(cowplot)

final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
output_name <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene/Spatial_BoxPlot_11Genes.pdf"
ggsave(output_name, plot = final_plot,
       width = 14, height = 4,
       units = c('in'), limitsize = FALSE,
       dpi = 300)


