##### Back to EdgeR!! 
#conda activate r4-base --> it has some issues in saving pdf files.. ## Seems like it's sapelo2 issue..
# conda activate r_env
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)
library("optparse")

option_list = list(
  make_option(c("--S1_Sparse"), type="character",
              help="S1_Sparse", metavar="character"),
  make_option(c("--S2_Sparse"), type="character",
              help="S2_Sparse"),
  make_option(c("--S1Name"), type="character",
              help="S1Name", metavar="character"),
  make_option(c("--S2Name"), type="character",
              help="S2Name"),
  make_option(c("--S1_Meta"), type="character",
              help="S1_Meta", metavar="character"),
  make_option(c("--S2_Meta"), type="character",
              help="S2_Meta", metavar="character"),
  make_option(c("--S1and2_500bpPeak"), type="character",
              help="S1and2_500bpPeak", metavar="character"),
  make_option(c("--S1and2_500bpInterPeak"), type="character",
              help="S1and2_500bpInterPeak", metavar="character"),
  make_option(c("--ClusterColumnName"), type="character",
              help="ClusterColumnName", metavar="character"),
  make_option(c("--OutFileName"), type="character",
              help="OutFileName", metavar="character"),
  make_option(c("--OutPath"), type="character",
              help="OutPath", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

S1_Name <- opt$S1Name
S2_Name <- opt$S2Name

Sparsefile_A619 <- opt$S1_Sparse
Sparsefile_Bif3 <- opt$S2_Sparse
MetaFileA619 <- opt$S1_Meta
MetaFileBif3 <- opt$S2_Meta

Peak_AllFile <- opt$S1and2_500bpPeak 
InterGenicFile <- opt$S1and2_500bpInterPeak

cluster_name <- opt$ClusterColumnName

OutfileName <- opt$OutFileName
WD <-   opt$OutPath
## 1. Load files and get filtered Sparse file for save space
## --> should combine A619+Bif3 in the begining to keep all the peaks as features

#Sparsefile_Bif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/Bif3_toComPeak.sparse"
#Sparsefile_A619  <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/A619_toComPeak.sparse"
#MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
#MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
Peak_All <-read.table(Peak_AllFile,header=F)
dim(Peak_All)
Peak_All_Pos <- paste(Peak_All$V1,Peak_All$V2,Peak_All$V3,sep="_")
head(Peak_All_Pos)
length(Peak_All_Pos)
#Sparse<- read.table(Sparsefile_A619,header=F)
#head(Sparse)

GetFilteredSparseData <- function(Sparsefile,MetaFile,
                                  SelectedPeaksPos,cluster_name){
  print("Start")
  #Sparsefile <- Sparsefile_Bif3
  #MetaFile <- MetaFileBif3
  #SelectedPeaksPos <- Peak_All_Pos
  Sparse<- read.table(Sparsefile,header=F)
  print("Done with reading data")
  meta_data <- read.table(MetaFile, header=TRUE)
  colnames(Sparse) <- c("PeakLocus","cellID","accessability")
  Sparse_Selected <- Sparse[Sparse$PeakLocus%in%SelectedPeaksPos,]
  Sparse_SelectedCells <-Sparse_Selected[Sparse_Selected$cellID %in% meta_data$cellID,]

  Sparse_SelectedCells$Celltype <- "Temp"
  celltypes <- levels(as.factor(meta_data[,cluster_name]))
  celltype <- celltypes[1]
  for (celltype in celltypes){
    print("Celltypes:")
    print(celltype)
    SelectedCellBarcode <- rownames(meta_data[which(meta_data[,cluster_name] == celltype),])
    Sparse_SelectedCells[Sparse_SelectedCells$cellID %in% SelectedCellBarcode,]$Celltype <- celltype
  }
  return(Sparse_SelectedCells)
}

setwd(WD)

Sparse_A619 <- GetFilteredSparseData(Sparsefile_A619,MetaFileA619,Peak_All_Pos,cluster_name)
Sparse_Bif3 <- GetFilteredSparseData(Sparsefile_Bif3,MetaFileBif3,Peak_All_Pos,cluster_name)
Sparse_A619$Celltype <- paste(S1_Name,Sparse_A619$Celltype,sep="_")
Sparse_Bif3$Celltype <- paste(S2_Name,Sparse_Bif3$Celltype,sep="_")

Sparse_Sum <- rbind(Sparse_A619,Sparse_Bif3)
head(Sparse_Sum)
dim(Sparse_Sum)

A619_Bif3_Celltype <- Sparse_Sum %>%
  group_by(Celltype,PeakLocus) %>%
  summarise_at(c("accessability"), sum)
head(A619_Bif3_Celltype)

A619_Bif3_Celltype_Count <-  spread(A619_Bif3_Celltype,key = Celltype,value =accessability)
tail(A619_Bif3_Celltype_Count)

saveRDS(A619_Bif3_Celltype_Count, file=paste0(OutfileName,"_AllPeaks_perCellType_Counts_FilteringBadPeaks.rds"))

## Note here!!!!

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/1.Correlation")
A619_Bif3_Celltype_Count <- readRDS("A619andrel2_AllPeaks_perCellType_Counts_FilteringBadPeaks.rds")
head(A619_Bif3_Celltype_Count)
InterGenicFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks_Intergenic.bed"
InterGenic <- read.table(InterGenicFile)
dim(InterGenic)
Intergenic_pos <- paste(InterGenic$V1,InterGenic$V2,InterGenic$V3,sep="_")
head(Intergenic_pos)
length(Intergenic_pos)
dim(A619_Bif3_Celltype_Count)
A619_Bif3_Celltype_Count_InterGenic <- A619_Bif3_Celltype_Count[A619_Bif3_Celltype_Count$PeakLocus%in%Intergenic_pos,]
dim(A619_Bif3_Celltype_Count_InterGenic)
## 2.Normalization
### Alt CPM Calc ## Normalization together with A619 and Bif3
PeakInfo <- A619_Bif3_Celltype_Count_InterGenic$PeakLocus
str(PeakInfo)
A619_Bif3_Celltype_Count <- A619_Bif3_Celltype_Count_InterGenic[,-1]
A619_Bif3_Celltype_Count[is.na(A619_Bif3_Celltype_Count)] <- 0

## 1) CPM normalization
head(A619_Bif3_Celltype_Count)
DevidedbySum <- apply(A619_Bif3_Celltype_Count,2,function(x){x/sum(x)})
A619_Bif3_CPM <- DevidedbySum*1000000

## 2) QuantileNormalization


head(A619_Bif3_CPM)
dim(A619_Bif3_CPM)

## This should be added!!!
count_A619 <- length(grep(paste0("^",S1_Name), colnames(A619_Bif3_CPM)))
A619_CPM <- A619_Bif3_CPM[,c(1:count_A619)]
Bif3_CPM <- A619_Bif3_CPM[,c((count_A619+1):ncol(A619_Bif3_CPM))]
S1Name_aligned <- gsub(S1_Name, S2_Name, colnames(A619_CPM))
get_order <- function(item, aligned_vector_A, original_vector_B) {
  if (item %in% aligned_vector_A) {
    return(match(item, aligned_vector_A))
  } else {
    return(max(length(aligned_vector_A), match(item, original_vector_B)))
  }
}
ordered_S2Name <- colnames(Bif3_CPM)[order(sapply(colnames(Bif3_CPM), get_order, S1Name_aligned, colnames(Bif3_CPM)))]
Bif3_CPM_Ordered <- Bif3_CPM[,ordered_S2Name]
head(Bif3_CPM_Ordered)
library(preprocessCore)

head(as.matrix(A619_Bif3_CPM))
A619_Bif3_Quantile <- normalize.quantiles(A619_Bif3_CPM)
head(A619_Bif3_Quantile)
dim(A619_Bif3_Quantile)
colnames(A619_Bif3_Quantile) <- colnames(A619_Bif3_CPM)
head(A619_Bif3_Quantile)

## This should be added!!!
A619_Q <- A619_Bif3_Quantile[,c(1:count_A619)]
Bif3_Q <- A619_Bif3_Quantile[,c((count_A619+1):ncol(A619_Bif3_CPM))]
dim(Bif3_Q)
Correlation <- cor(A619_Q,Bif3_Q,  method = "pearson")

## 3) Get the most variable 2000 ACR.
head(A619_Bif3_Quantile)
Variance_by_Peak <- apply(A619_Bif3_Quantile,1,var)
length(Variance_by_Peak)
length(Intergenic_pos)
sort(Variance_by_Peak, decreasing = TRUE)[2000]
A619_Bif3_Quantile_Top2000 <- A619_Bif3_Quantile[Variance_by_Peak >= sort(Variance_by_Peak, decreasing = TRUE)[2000],]
dim(A619_Bif3_Quantile_Top2000)
head(A619_Bif3_Quantile_Top2000)
## This should be added!!!
A619_Q_Top2000 <- A619_Bif3_Quantile_Top2000[,c(1:count_A619)]
Bif3_Q_Top2000 <- A619_Bif3_Quantile_Top2000[,c((count_A619+1):ncol(A619_Bif3_CPM))]
colnames(A619_Q_Top2000)
colnames(Bif3_Q_Top2000)
Correlation_Top2000 <- cor(A619_Q_Top2000,Bif3_Q_Top2000,  method = "pearson")
Correlation_Top2000

### 4) Visulaize the plot
library(reshape2)
library(ggplot2)

melted_Correlation <- melt(Correlation, na.rm = TRUE)
head(melted_Correlation)

ggplot(data = melted_Correlation, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "#2de361", high = "red", mid =  "#e3e32d",
                       midpoint = mean(melted_Correlation$value), 
                       limit = c(min(melted_Correlation$value),max(melted_Correlation$value)), space = "Lab", 
                       name="Pearson\nCorrelation")+
  
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() +  xlab(S1_Name) + ylab(S2_Name)

ggsave(paste0(OutfileName,"_IntergenicPeaks.pdf"), width=10, height=10)

###
#pdf(file = NULL)
melted_Correlation_Top2000 <- melt(Correlation_Top2000, na.rm = TRUE)
head(melted_Correlation_Top2000)

ggplot(data = melted_Correlation_Top2000, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "#002aff", high = "#b30404", mid =  "white",
                       midpoint = median(melted_Correlation_Top2000$value), 
                       limit = c(min(melted_Correlation_Top2000$value),max(melted_Correlation_Top2000$value)), space = "Lab", 
                       name="Pearson\nCorrelation")+
  #scale_fill_gradient(low = "#002aff", high = "red",
  #                     limit = c(min(melted_Correlation_Top2000$value),max(melted_Correlation_Top2000$value)), space = "Lab", 
  #                     name="Pearson\nCorrelation")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() +  xlab(S1_Name) + ylab(S2_Name)

ggsave(paste0(OutfileName,"_IntergenicPeaks_Top2000ACR.pdf"), width=10, height=10)

###
#pdf(file = NULL)
Correlation_Top2000_A619 <- cor(A619_Q_Top2000,A619_Q_Top2000,  method = "pearson")

melted_Correlation_Top2000_A619 <- melt(Correlation_Top2000_A619, na.rm = TRUE)
head(melted_Correlation_Top2000_A619)

ggplot(data = melted_Correlation_Top2000_A619, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "#002aff", high = "#b30404", mid =  "white",
                       midpoint = median(melted_Correlation_Top2000_A619$value), 
                       limit = c(min(melted_Correlation_Top2000_A619$value),max(melted_Correlation_Top2000_A619$value)), space = "Lab", 
                       name="Pearson\nCorrelation")+
  #scale_fill_gradient(low = "#002aff", high = "red",
  #                     limit = c(min(melted_Correlation_Top2000$value),max(melted_Correlation_Top2000$value)), space = "Lab", 
  #                     name="Pearson\nCorrelation")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() +  xlab(S1_Name) + ylab(S2_Name)

ggsave(paste0(OutfileName,"_IntergenicPeaks_Top2000ACR.pdf"), width=10, height=10)

