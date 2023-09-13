##### Back to EdgeR!! 
#conda activate r4-base --> it has some issues in saving pdf files.. ## Seems like it's sapelo2 issue..
# conda activate r_env
##### I will calculate ACRs only in intergenic regions.
### By Celltype.

library(edgeR)
library(tidyverse)
library(stringr)
library("optparse")
library(rlang)

option_list = list(
  make_option(c("--WD"), type="character", 
              help="WD", metavar="character"),
  make_option(c("--Sparse"), type="character", 
              help="Sparse", metavar="character"),
  make_option(c("--Meta"), type="character",
              help="Meta", metavar="character"),
  make_option(c("--PeakI"), type="character", 
              help="PeakIntergenic", metavar="character"),
  make_option(c("--PeakG"), type="character", 
              help="PeakGenic", metavar="character"),
  make_option(c("--SampleName"), type="character", 
              help="SampleName", metavar="character")
  
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Sparsefile_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks.sparse"
#MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
#SampleName <- "A619"
#Peak_Inter <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic.bed",header=F)
#Peak_Genic <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Genic.bed",header=F)

#Sparsefile_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/bif3.500bp_peaks.sparse"
#MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
#SampleName <- "bif3"
#Peak_Inter <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/bif3.500bp_peaks_Intergenic.bed",header=F)
#Peak_Genic <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/bif3.500bp_peaks_Genic.bed",header=F)


Sparsefile_A619 <- opt$Sparse
MetaFileA619 <- opt$Meta
SampleName <- opt$SampleName
Peak_Inter <- read.table(opt$PeakI,header=F)
Peak_Genic <- read.table(opt$PeakG,header=F)
WD <- opt$WD
## 1. Load files and get filtered Sparse file for save space
## --> should combine A619+Bif3 in the begining to keep all the peaks as features

dim(Peak_Genic)+dim(Peak_Inter)
Peak_All <- rbind(Peak_Genic,Peak_Inter)
head(Peak_All)

Peak_All_Pos <- paste(Peak_All$V1,Peak_All$V2,Peak_All$V3,sep="_")
head(Peak_All_Pos)
length(Peak_All_Pos)
print(Sparsefile_A619)
#Sparse<- read.table(Sparsefile_A619,header=F)
#head(Sparse)
#length(levels(as.factor(Sparse$V1)))

GetFilteredSparseData <- function(Sparsefile,MetaFile,
                                  SelectedPeaksPos,cluster_name="Ann_v3"){
  print("Start")
  #Sparsefile <- Sparsefile_A619
  #MetaFile <- MetaFileA619
  #SelectedPeaksPos <- Peak_All_Pos
  Sparse<- read.table(Sparsefile,header=F)
  print("Done with reading data")
  meta_data <- read.table(MetaFile, header=TRUE)
  #print(head(meta_data))
  colnames(Sparse) <- c("PeakLocus","cellID","accessability")
  Sparse_Selected <- Sparse[Sparse$PeakLocus%in%SelectedPeaksPos,]
  #print(head(Sparse_Selected))
  Sparse_SelectedCells <-Sparse_Selected[Sparse_Selected$cellID %in% meta_data$cellID,]
  #print(head(Sparse_SelectedCells))
  head(meta_data)
  Sparse_SelectedCells$Celltype <- "Temp"
  celltypes <- levels(as.factor(meta_data[,cluster_name]))
  celltype <- celltypes[1]
  for (celltype in celltypes){
    print(celltype)
    SelectedCellBarcode <- rownames(meta_data[which(meta_data[,cluster_name] == celltype),])
    Sparse_SelectedCells[Sparse_SelectedCells$cellID %in% SelectedCellBarcode,]$Celltype <- celltype
  }
  return(Sparse_SelectedCells)
}
print(MetaFileA619)
Sparse_A619 <- GetFilteredSparseData(Sparsefile_A619,MetaFileA619,Peak_All_Pos)
head(Sparse_A619)

Sparse_A619[endsWith(Sparse_A619$cellID, '_2'),]$Celltype <- paste0("Re2_",
                    Sparse_A619[endsWith(Sparse_A619$cellID, '_2'),]$Celltype)
Sparse_A619[!endsWith(Sparse_A619$cellID, '_2'),]$Celltype <- paste0("Re1_",
                    Sparse_A619[!endsWith(Sparse_A619$cellID, '_2'),]$Celltype)

tail(Sparse_A619)

A619_Celltype <- Sparse_A619 %>%
  group_by(Celltype,PeakLocus) %>%
  summarise_at(c("accessability"), sum)
head(A619_Celltype)

A619_Celltype_Count <-  spread(A619_Celltype,key = Celltype,value =accessability)
tail(A619_Celltype_Count)

setwd(WD)
#setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/9.Correlation")
saveRDS(A619_Celltype_Count, file=paste0(SampleName,"_Replicates.rds"))

## Note here!!!!
#A619_Celltype_Count <- readRDS("bif3_Replicates.rds")
head(A619_Celltype_Count)

#A619_Bif3_Celltype_Count_InterGenic <- A619_Bif3_Celltype_Count[A619_Bif3_Celltype_Count$PeakLocus%in%Intergenic_pos,]
#dim(A619_Bif3_Celltype_Count_InterGenic)

## 2.Normalization
### Alt CPM Calc
PeakInfo <- A619_Celltype_Count$PeakLocus
str(PeakInfo)
A619_Celltype_Count <- A619_Celltype_Count[,-1]
A619_Celltype_Count[is.na(A619_Celltype_Count)] <- 0

## 1) CPM normalization
head(A619_Celltype_Count)
DevidedbySum <- apply(A619_Celltype_Count,2,function(x){x/sum(x)})
A619_CPM <- DevidedbySum*1000000

## 2) QuantileNormalization


head(A619_CPM)
dim(A619_CPM)
colnames(A619_CPM)
A619_Re1 <- A619_CPM[,c(1:13)]
A619_Re2 <- A619_CPM[,c(14:26)]

library(preprocessCore)
head(as.matrix(A619_CPM))
A619_Quantile <- normalize.quantiles(A619_CPM)
head(A619_Quantile)
colnames(A619_Quantile) <- colnames(A619_CPM)
head(A619_Quantile)
A619_Re1_Q <- A619_Quantile[,c(1:13)]
A619_Re2_Q <- A619_Quantile[,c(14:26)]
dim(A619_Re1_Q)
Correlation <- cor(A619_Re1_Q,A619_Re2_Q,  method = "pearson")

## 3) Get the most variable 2000 ACR.
head(A619_Quantile)
rownames(A619_Quantile) <- PeakInfo
Variance_by_Peak <- apply(A619_Quantile,1,var)
length(Variance_by_Peak)
sort(Variance_by_Peak, decreasing = TRUE)[2000]
A619_Quantile_Top2000 <- A619_Quantile[Variance_by_Peak >= sort(Variance_by_Peak, decreasing = TRUE)[2000],]
dim(A619_Quantile_Top2000)
head(A619_Quantile_Top2000)
Top2000_Row <- rownames(A619_Quantile_Top2000)
head(A619_Celltype_Count)
# Convert tibble to data.frame
A619_Celltype_Count_df <- as.data.frame(A619_Celltype_Count)
# Set the row names
rownames(A619_Celltype_Count_df) <- PeakInfo
#A619_Celltype_Count_df[Top2000_Row,][,6]
#A619_Celltype_Count_df[Top2000_Row,][,7]
#head(A619_Celltype_Count_df)
sum(A619_Celltype_Count_df[Top2000_Row,][,6]==0)
sum(A619_Celltype_Count_df[Top2000_Row,][,7]==0)
sum(A619_Celltype_Count_df[Top2000_Row,][,8]==0)
head(A619_Celltype_Count_df)
dim(A619_Celltype_Count_df)
A619_Celltype_Count_df_Filtered <- A619_Celltype_Count_df[Top2000_Row,][A619_Celltype_Count_df[Top2000_Row,][,6]!=0,]
dim(A619_Celltype_Count_df_Filtered)

A619_Q_Top2000_Re1 <- A619_Quantile_Top2000[,c(1:13)]
A619_Q_Top2000_Re2 <- A619_Quantile_Top2000[,c(14:26)]
Correlation_Top2000 <- cor(A619_Q_Top2000_Re1,A619_Q_Top2000_Re2,  method = "pearson")
#Correlation_Top2000

### 4) Visulaize the plot
library(reshape2)
library(ggplot2)


CorrPlot_Function <- function(Correlation,SampleName,Prefix=".pdf"){
  melted_Correlation <- melt(Correlation, na.rm = TRUE)
  head(melted_Correlation)
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
  coord_fixed() +  xlab("Replicate1") + ylab("Replicate2")

ggsave(paste0(SampleName,Prefix), width=10, height=10)
}

CorrPlot_Function(Correlation,SampleName,Prefix="_Replicates_AllACRs.pdf")
CorrPlot_Function(Correlation_Top2000,SampleName,Prefix="_Replicates_Top2000Variable.pdf")
