library(edgeR)
library(ggplot2)

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


OCCells <- read.table("/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/DEGTestResult_MannuallySelectedCells.csv",
                      sep = ",", header = TRUE,  quote = "")
TargetGeneFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/dACR_withTAATInfo.txt"
TargetGene <- read.table(TargetGeneFile,fill=TRUE,header=TRUE)
TargetGene_Filtered <- TargetGene[TargetGene$TAAT =="TAAT",]
TargeneGeneList <- TargetGene_Filtered$gene_model

TargeneGeneList
length(TargeneGeneList)
NonTargetGeneList <- setdiff(Result_table$genes, TargeneGeneList)
RandomGenes <- sample(NonTargetGeneList, size = 381, replace = FALSE)
Genes_excluded <- setdiff(A, B)
TargetGenes <- OCCells[OCCells$genes %in% TargeneGeneList, ]$logFC
NontargetGenes <- OCCells[OCCells$genes %in% RandomGenes, ]$logFC
mean(TargetGenes)
mean(NontargetGenes)
median(TargetGenes)
median(NontargetGenes)
gene_data <- data.frame(
  Expression = c(TargetGenes, NontargetGenes),
  GeneType = c(rep("TAATMotifGenes", length(TargetGenes)), rep("RandomGenes", length(NontargetGenes)))
)

# Create the box plot using ggplot2
ggplot(gene_data, aes(x = GeneType, y = Expression, fill = GeneType)) +
  geom_violin(trim = FALSE) +
  labs(title = " ", y = "logFC") +
  theme_minimal()


ggsave("/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/TAATmotif_violetplot_spatialRNAseq.pdf", width = 8, height = 6)

