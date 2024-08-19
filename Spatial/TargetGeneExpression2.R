library(edgeR)
library(ggplot2)

file_path <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene/WTFive_Bif3Ten_count.csv"
data <- read.csv(file_path, row.names = 1, header = TRUE)
#data <- as.data.table(data)
dim(data)
group <- factor(c(rep("WT",5), rep("Bif3",10)))
dge <- DGEList(counts=data, group=group)
dge <- calcNormFactors(dge, method="TMM")
tmm_normalized_counts <- cpm(dge, normalized.lib.sizes=TRUE)
head(tmm_normalized_counts)

AllPlotgenes_symbol <- c('arftf4', 'arftf30', 'arftf18', 'arftf3', 'arftf20', 'arftf25', 'arftf10', 'arftf36', 'arftf23', 'arftf26', 'knox1')
AllPlotgenes <- c('Zm00001eb433460', 'Zm00001eb066640', 'Zm00001eb067270', 'Zm00001eb142540', 'Zm00001eb224680', 'Zm00001eb232120', 'Zm00001eb243930', 'Zm00001eb292830', 'Zm00001eb363810', 'Zm00001eb370810', 'Zm00001eb001720')

tmm_normalized_selected <- tmm_normalized_counts[AllPlotgenes,]

Plotlist <- list()
i <-11
for (i in c(1:nrow(tmm_normalized_selected))) {
  # Extract gene data for WT and Bif3
  gene <- rownames(tmm_normalized_selected)[[i]]
  gene_data <- tmm_normalized_selected[gene, ]
  wt_values <- gene_data[c(1:5)]
  bif3_values <- gene_data[c(6:15)]
  
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
output_name <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene/VioletPlot_11Genes.pdf"
ggsave(output_name, plot = final_plot,
       width = 14, height = 4,
       units = c('in'), limitsize = FALSE,
       dpi = 300)
