library(edgeR)
library(ggplot2)
library(cowplot)

##### 1) load Data and write the Replicate info!
file_path <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/"

# WT Re1
data <- read.csv(paste0(file_path,"A619_B_MannuallySelected.csv"), row.names = 1, header = TRUE)
#data <- as.data.table(data)
WT_Re1 <- colSums(data)
#WT Re2
data <- read.csv(paste0(file_path,"A619_C_MannuallySelected.csv"), row.names = 1, header = TRUE)
WT_Re2_1 <- colSums(data)
data <- read.csv(paste0(file_path,"A619_D_WT_MannuallySelected.csv"), row.names = 1, header = TRUE)
WT_Re2_2 <- colSums(data)
WT_Re2 <- WT_Re2_1+WT_Re2_2

# Bif3 Re1
data <- read.csv(paste0(file_path,"A619_D_Bif3_Cluser6.csv"), row.names = 1, header = TRUE)
Bif3_Re1 <- colSums(data)

# Bif3 Re2
data <- read.csv(paste0(file_path,"Bif3_A_Cluser1.csv"), row.names = 1, header = TRUE)
Bif3_Re2 <- colSums(data)

# Bif3 Re3

data <- read.csv(paste0(file_path,"Bif3_B_Cluser6.csv"), row.names = 1, header = TRUE)
Bif3_Re3 <- colSums(data)

# Bif3 Re4
data <- read.csv(paste0(file_path,"Bif3_C_Cluser11.csv"), row.names = 1, header = TRUE)
Bif3_Re4_1 <- colSums(data)
data <- read.csv(paste0(file_path,"Bif3_D_Cluser7.csv"), row.names = 1, header = TRUE)
Bif3_Re4_4 <- colSums(data)
Bif3_Re4 <- Bif3_Re4_1+Bif3_Re4_4

CountTable <- data.frame(WT_Re1,WT_Re2,Bif3_Re1,Bif3_Re2,Bif3_Re3,Bif3_Re4)
CountTable <- CountTable[grep("^Zm00", rownames(CountTable)), ]
head(CountTable)
tail(CountTable)

### EdgeR count info: row - gene column Sample
group <- factor(c(rep("WT",2), rep("Bif3",4)))
SampleID <- colnames(CountTable)
targets <- data.frame(SampleID ,group)
targets$group <- factor(targets$group)
targets$group <- factor(targets$group,levels=c("WT","Bif3"))
design <- model.matrix(~group, data=targets)
log_cpm <- cpm(CountTable, log = TRUE)
y <- DGEList(counts=CountTable, gene=rownames(CountTable))
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMRobustDisp(y,design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit,coef = 2)
Result_table <- topTags(lrt, n=dim(CountTable)[1], sort.by="none")$table
head(Result_table)

tmm_normalized_counts <- cpm(y, normalized.lib.sizes=TRUE,log=T)
write.table(Result_table, file = "/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/DEGTestResult_MannuallySelectedCells.csv",
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


AllPlotgenes_symbol <- c('arftf4', 'arftf23','arftf30', 'arftf18',
                         'arftf3', 'arftf20', 'arftf25',
                         'arftf10', 'arftf36', 
                         'arftf26', 'knox1')
AllPlotgenes <- c('Zm00001eb067270','Zm00001eb292830', 'Zm00001eb433460', 'Zm00001eb224680',
                  'Zm00001eb066640', 'Zm00001eb232120', 'Zm00001eb363810', 
                  'Zm00001eb142540', 'Zm00001eb243930',  
                  'Zm00001eb370810', 'Zm00001eb001720')

tmm_normalized_selected <- tmm_normalized_counts[AllPlotgenes,]
RawCount_Selected <- CountTable[AllPlotgenes,]
RawCount_Selected
head(tmm_normalized_selected)
Plotlist <- list()
i <-1
i <- 9
head(log_cpm)
log_cpm_selected <- log_cpm[AllPlotgenes,]
AllPlotgenes_symbol
log_cpm_selected2 <- data.frame(log_cpm_selected)
log_cpm_selected2$geneName <- c(AllPlotgenes_symbol)

write.table(log_cpm_selected2, file = "/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/CPMValues_TargetGene.csv",
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


for (i in c(1:nrow(log_cpm_selected))) {
  # Extract gene data for WT and Bif3
  gene <- rownames(log_cpm_selected)[[i]]
  gene_data <- log_cpm_selected[gene, ]
  wt_values <- gene_data[c(1:2)]
  bif3_values <- gene_data[c(3:6)]
  FDR <- round(Result_table[gene,]$FDR,5)
  PValue <- round(Result_table[gene,]$PValue,5)
  # Create a data frame for ggplot
  plot_data <- data.frame(
    value = c(wt_values, bif3_values),
    group = rep(c("WT", "Bif3"), c(length(wt_values), length(bif3_values)))
  )
  plot_data$group <- factor(plot_data$group,levels=c("WT","Bif3"))
  # Plotting the violin plot
  p <- ggplot(plot_data, aes(x = group, y = value)) +
    geom_boxplot() +
    labs(title = paste0(AllPlotgenes_symbol[[i]],"\nPvalue:",PValue,"\nFDR:",FDR),
         x = "Group", y = "log CPM") +
    theme_minimal()+
    ylim(0, NA)
  Plotlist[[i]] <- p
}

library(cowplot)

final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
output_name <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/LogCPM_Spatial_BoxPlot_11Genes_MannuallySelectedCells.pdf"
ggsave(output_name, plot = final_plot,
       width = 14, height = 4,
       units = c('in'), limitsize = FALSE,
       dpi = 300)




