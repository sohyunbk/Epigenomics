##### Back to EdgeR!!
#conda activate r4-base --> it has some issues in saving pdf files.. ## Seems like it's sapelo2 issue..
# conda activate r_env
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)
library("optparse")
library(dplyr) 
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
              help="OutPath", metavar="character"),
  make_option(c("--CellTypeOrder"), type="character",
          help="CellTypeOrder", metavar="character")
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
CTOrders <-   opt$CellTypeOrder

CTOrder <- readLines(CTOrders)

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
Sparse_A619$Celltype <- paste0(S1_Name,Sparse_A619$Celltype)
Sparse_Bif3$Celltype <- paste0(S2_Name,Sparse_Bif3$Celltype)

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
#A619_Bif3_Celltype_Count <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/1.Correlation/A619andBif3_CTNameReverse_AllPeaks_perCellType_Counts_FilteringBadPeaks.rds")
#InterGenicFile <-"/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4//A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks_Intergenic.bed"
## Note here!!!!
#A619_Bif3_Celltype_Count <- readRDS("AllPeaks_perCellType_Counts.rds")
head(A619_Bif3_Celltype_Count)


InterGenic <- read.table(InterGenicFile)
dim(InterGenic)
Intergenic_pos <- paste(InterGenic$V1,InterGenic$V2,InterGenic$V3,sep="_")
head(Intergenic_pos)
length(Intergenic_pos)
dim(A619_Bif3_Celltype_Count)
A619_Bif3_Celltype_Count_InterGenic <- A619_Bif3_Celltype_Count[A619_Bif3_Celltype_Count$PeakLocus%in%Intergenic_pos,]
dim(A619_Bif3_Celltype_Count_InterGenic)
## 2.Normalization
### Alt CPM Calc
PeakInfo <- A619_Bif3_Celltype_Count_InterGenic$PeakLocus
str(PeakInfo)
A619_Bif3_Celltype_Count <- A619_Bif3_Celltype_Count_InterGenic[,-1]
A619_Bif3_Celltype_Count[is.na(A619_Bif3_Celltype_Count)] <- 0

## 1) CPM normalization
#head(A619_Bif3_Celltype_Count)
DevidedbySum <- apply(A619_Bif3_Celltype_Count,2,function(x){x/sum(x)})
A619_Bif3_CPM <- DevidedbySum*1000000

## 2) QuantileNormalization
library(preprocessCore)
#head(as.matrix(A619_Bif3_CPM))
A619_Bif3_Quantile <- normalize.quantiles(A619_Bif3_CPM)
#head(A619_Bif3_Quantile)
dim(A619_Bif3_Quantile)
colnames(A619_Bif3_Quantile) <- colnames(A619_Bif3_CPM)
#head(A619_Bif3_Quantile)

count_A619 <- length(grep(paste0("^",S1_Name), colnames(A619_Bif3_CPM)))
A619_Q <- A619_Bif3_Quantile[,c(1:count_A619)]
Bif3_Q <- A619_Bif3_Quantile[,c((count_A619+1):ncol(A619_Bif3_CPM))]
S1Name_aligned <- gsub(S1_Name, S2_Name, colnames(A619_Q))

#CTOrder
print("Bif3_Q object")
print(CTOrder)
head(Bif3_Q)
Bif3_Q <- Bif3_Q[,paste0("relk1",CTOrder)]
A619_Q <- A619_Q[,paste0("A619",CTOrder)]
Correlation <- cor(A619_Q,Bif3_Q,  method = "pearson")

## 3) Get the most variable 2000 ACR.
#head(A619_Bif3_Quantile)
Variance_by_Peak <- apply(A619_Bif3_Quantile,1,var)
length(Variance_by_Peak)
length(Intergenic_pos)
sort(Variance_by_Peak, decreasing = TRUE)[2000]
A619_Bif3_Quantile_Top2000 <- A619_Bif3_Quantile[Variance_by_Peak >= sort(Variance_by_Peak, decreasing = TRUE)[2000],]
dim(A619_Bif3_Quantile_Top2000)
#head(A619_Bif3_Quantile_Top2000)

A619_Q_Top2000 <- A619_Bif3_Quantile_Top2000[,c(1:count_A619)]
Bif3_Q_Top2000 <- A619_Bif3_Quantile_Top2000[,c((count_A619+1):ncol(A619_Bif3_CPM))]
## I should edit this part!
Bif3_Q_Top2000 <- Bif3_Q_Top2000[,paste0("relk1",CTOrder)]
A619_Q_Top2000 <- A619_Q_Top2000[,paste0("A619",CTOrder)]

cols_to_keep <- !grepl("Unknown|G2", colnames(Bif3_Q_Top2000))
Bif3_Q_Top2000 <- Bif3_Q_Top2000[, cols_to_keep]
A619_Q_Top2000 <- A619_Q_Top2000[, cols_to_keep]


Correlation_Top2000 <- cor(A619_Q_Top2000,Bif3_Q_Top2000,  method = "pearson")
#Correlation_Top2000 <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/1.Correlation/A619andBif3_CTNameReverse_Top2000Correlation.txt")
write.table(Correlation_Top2000,
            file = paste0(OutfileName,"_Top2000Correlation.txt"), sep = "\t", 
            row.names = TRUE, quote=F, col.names = TRUE)
### 4) Visulaize the plot
library(reshape2)
library(ggplot2)
######### Draw Plot
CorrPlot_Function <- function(Correlation,SampleName,Prefix=".pdf"){
  melted_Correlation <- melt(Correlation, na.rm = TRUE)

  #head(melted_Correlation)
  Min <- min(melted_Correlation$value)
  Max <- max(melted_Correlation$value)

  ggplot(data = melted_Correlation, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()+
    scale_fill_gradient2(low = "#002aff", high = "red", mid =  "white",
                         midpoint = (Min+Max)/2,
                         limit = c(Min,Max),
                         breaks=c(round(Min,2)+0.01,round((Min+Max)/2,2),round(Max,2)-0.01),
                         name="Pearson\nCorrelation")+

    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1))+
    coord_fixed() +  xlab(S1_Name) + ylab(S2_Name)

  ggsave(paste0(SampleName,Prefix), width=10, height=10)
}
CorrPlot_Function(Correlation,OutfileName,Prefix="_AllACRs.pdf")
CorrPlot_Function(Correlation_Top2000,OutfileName,Prefix="_Top2000Variable.pdf")
