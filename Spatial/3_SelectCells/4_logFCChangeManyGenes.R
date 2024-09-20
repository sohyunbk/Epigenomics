library(edgeR)
library(ggplot2)


OCCells <- read.table("/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/DEGTestResult_MannuallySelectedCells.csv",
                      sep = ",", header = TRUE,  quote = "")
TargetGeneFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/dACR_withTAATInfo.txt"
TargetGene <- read.table(TargetGeneFile,fill=TRUE,header=TRUE)
TargetGene_Filtered <- TargetGene[TargetGene$TAAT =="TAAT",]
TargeneGeneList <- TargetGene_Filtered$gene_model

TargeneGeneList
length(TargeneGeneList)
NonTargetGeneList <- setdiff(OCCells$genes, TargeneGeneList)
RandomGenes <- sample(NonTargetGeneList, size = 381, replace = FALSE)
TargetGenes <- OCCells[OCCells$genes %in% TargeneGeneList, ]$logFC
NontargetGenes <- OCCells[OCCells$genes %in% RandomGenes, ]$logFC
mean(TargetGenes)
mean(NontargetGenes)
NontargetGenes$GeneType <- 
gene_data <- data.frame(
  Expression = c(TargetGenes, NontargetGenes),
  GeneType = c(rep("TAATMotifGenes", length(TargetGenes)), rep("RandomGenes", length(NontargetGenes)))
)

# Create the box plot using ggplot2
PlotList <- list()
PlotList[[1]] <- ggplot(gene_data, aes(x = GeneType, y = Expression, group = GeneType)) +
  geom_boxplot(fill=NA) +
  labs(title = " ", y = "logFC") +
  theme_minimal()

#ggsave("/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/TAATmotif_violetplot_spatialRNAseq.pdf", width = 8, height = 6)

######## Decreased gene expression ##################################
Dir= "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC/"

MakePlot <- function(FileName,ConditionName){
TargetGeneFile <- paste0(Dir,FileName)
TargetGene <- read.table(TargetGeneFile,fill=TRUE,header=TRUE)
TargeneGeneList <-rownames(TargetGene)

NonTargetGeneList <- setdiff(OCCells$genes, TargeneGeneList)
RandomGenes <- sample(NonTargetGeneList, size = length(TargeneGeneList), replace = FALSE)
TargetGenes <- OCCells[OCCells$genes %in% TargeneGeneList, ]$logFC
NontargetGenes <- OCCells[OCCells$genes %in% RandomGenes, ]$logFC
print(ConditionName)
print(mean(TargetGenes))
print(mean(NontargetGenes))
gene_data <- data.frame(
  Expression = c(TargetGenes, NontargetGenes),
  GeneType = c(rep(ConditionName, length(TargetGenes)), rep("RandomGenes", length(NontargetGenes)))
)

plot <- ggplot(gene_data, aes(x = GeneType, y = Expression, group = GeneType)) +
  geom_boxplot(fill=NA) +
  labs(title = " ", y = "logFC") +
  theme_minimal()

return(plot)
}

PlotList[[2]] <- MakePlot("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_Fimo1_GCACAGCAGCR.txt","Fimo1_GCACAGCAGCR")
PlotList[[3]] <- MakePlot("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_Fimo3_GCAGCATGC.txt","Fimo3_GCAGCATGC")
PlotList[[4]] <- MakePlot("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_Fimo4_CGCGCCGCGCC.txt","Fimo4_CGCGCCGCGCC")
PlotList[[5]] <- MakePlot("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_Fimo5_YAGAGAGAGA.txt","Fimo5_YAGAGAGAGA")
PlotList[[6]] <- MakePlot("BoxPlot_GeneBodyAcc_ClosestGene_withdACRA619Higher_Fimo7_GCTAGCTAGC.txt","Fimo7_GCTAGCTAGC")

library(cowplot)

final_plot <- plot_grid(plotlist = PlotList, ncol = 6)
output_name <- "/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/LogFC_fordACRwithTAATforOtherFimos.pdf"
ggsave(output_name, plot = final_plot,
       width = 14, height = 4,
       units = c('in'), limitsize = FALSE,
       dpi = 300)
