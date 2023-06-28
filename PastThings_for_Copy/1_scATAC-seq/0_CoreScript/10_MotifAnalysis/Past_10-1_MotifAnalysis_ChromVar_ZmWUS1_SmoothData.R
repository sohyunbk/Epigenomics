#https://greenleaflab.github.io/chromVAR/articles/Introduction.html#example-data
#conda activate JASPAR_act
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
set.seed(2017)
library(preprocessCore)
library(tidyverse)
library(stringr)

library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
#library(varistran)
library(edgeR)
library(parallel)
library(png)
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)
library(RANN)

#### bring the sparse data and then make it to matrix ####
# Column should be cell type and row should be peak
## Nead Granges about peak regions

#marker_correlation <- function(MetaFile, 
#                               SparseData, 
#                               cluster_name = "Ann_v3"){
#  Sparse_Bif3 <- read.table(,header=F)
#  
#}

SampleName <- "A619"
SampleName <- "Bif3"

if (SampleName == "A619"){
  Sparse_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/A619_toComPeak.sparse"
  MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
  pcs <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt")
  meta_data <- read.table(MetaFileA619, header=TRUE)
  Sparse <- read.table(Sparse_A619,header=F)
  } else if (SampleName == "Bif3") {
    Sparse_Bif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/Bif3_toComPeak.sparse"
    MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
    pcs <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt")    
    Sparse <- read.table(Sparse_Bif3,header=F)  
    meta_data <- read.table(MetaFileBif3, header=TRUE)
   }

############################
head(meta_data)
colnames(Sparse) <- c("PeakLocus","cellID","accessability")


InterGenic_peak <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/ComA619Bif3.unique500bpPeaks_BLRemove_Intergenic.bed",header=F)
dim(InterGenic_peak)

#Genic_peak_Pos <- paste(Genic_peak$V1,Genic_peak$V2,Genic_peak$V3,sep="_")
#head(Genic_peak_Pos)
InterGenic_peak_Pos <- paste(InterGenic_peak$V1,InterGenic_peak$V2,InterGenic_peak$V3,sep="_")
head(InterGenic_peak_Pos)
Sparse_Intergenic <- Sparse[Sparse$PeakLocus%in%InterGenic_peak_Pos,]
head(Sparse_Intergenic)
dim(Sparse_Intergenic)
Sparse_Intergenic_SelectedCells <-Sparse_Intergenic[Sparse_Intergenic$cellID %in% meta_data$cellID,]
dim(Sparse_Intergenic_SelectedCells)

cluster_name = "Ann_v3"
## 1) Make count matrix for input
#Sparse <- as_tibble(Sparse_Intergenic_SelectedCells)
## Add here to split the peaks into intergenic & genic regions
Peak_Cell_Count <-  spread(Sparse_Intergenic_SelectedCells,key = cellID,value =accessability)
Peak_Cell_Count[is.na(Peak_Cell_Count)] <- 0
dim(Peak_Cell_Count)

## 2) Make Pos Grange Input
#https://www.geeksforgeeks.org/how-to-split-column-into-multiple-columns-in-r-dataframe/
PeakPos <- data.frame(Pos = Peak_Cell_Count$PeakLocus)
PeakPos[c('Chr','Start', 'End')] <- str_split_fixed(PeakPos$Pos, '_', 3)
str(PeakPos)
head(PeakPos)
dim(PeakPos)
peaks <-  GRanges(seqnames = PeakPos$Chr,
                       ranges = IRanges(start = as.numeric(PeakPos$Start),
                        end = as.numeric(PeakPos$End)))

## 3) Make ChromVar input 
#Peak_Cell_Count <- readRDS("Bif3_Peak_perCell_Counts.rds")
str(Peak_Cell_Count)
Peak_Cell_Count_Input <- as.matrix(Peak_Cell_Count[,-1])
str(Peak_Cell_Count_Input)
head(Peak_Cell_Count_Input)[,c(1:10)]
dim(Peak_Cell_Count_Input)

ChromVarInput_counts <- SummarizedExperiment(assays = 
                                          list(counts = Peak_Cell_Count_Input),
                                        rowRanges = peaks)

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/WithGoodPeak")
SampleName

saveRDS(Peak_Cell_Count, file=paste0(SampleName,"_Peak_perCell_Counts.rds"))
saveRDS(ChromVarInput_counts, file=paste0(SampleName,"_ChromVarInput.rds"))

### 4) Get motifs!! ## getMatrixSet from TFBSTools
#https://bioconductor.org/packages/devel/bioc/manuals/TFBSTools/man/TFBSTools.pdf #Page10
#conda activate MotifR !!
## Should use this file from GEM: HB67_WUS1_B73v5_Q30_default_finalBl_2.all.PFM_MEME.txt : TGAA motifs until m0, m1, m2
## It is corresponding to file:///Users/sohyun/Documents/2.SingleCellATAC/6.WUS_DAP-seqData/HB67_WUS1_B73v5_Q30_default_finalBl_LowestThreshold_Default/HB67_WUS1_B73v5_Q30_default_finalBl_outputs/HB67_WUS1_B73v5_Q30_default_finalBl_2.results.htm

#ChromVarInput_counts <- readRDS(paste0(SampleName,"_ChromVarInput.rds"))
#Peak_Cell_Count <- readRDS(paste0(SampleName,"_Peak_perCell_Counts.rds"))

str(ChromVarInput_counts)
## Remove!! chr6	181356755	181357256	1.10868807050461


library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(BSgenome)

## Get motif 
opts <- list()
#Zea mays Taxonomy ID: 4577 , Ara  Taxonomy ID: 3702
opts[["species"]] <- 3702 #NCBI tax IDs (9606).
opts[["all_versions"]] <- TRUE
PFMatrixList_Ara <- getMatrixSet(JASPAR2022, opts)
PFMatrixList_Ara[['MA']]
opts <- list()
#Zea mays Taxonomy ID: 4577 , Ara  Taxonomy ID: 3702
opts[["species"]] <- 4577 #NCBI tax IDs (9606).
opts[["all_versions"]] <- TRUE
PFMatrixList_Maize <- getMatrixSet(JASPAR2022, opts)
PFMatrixList <- c(PFMatrixList_Ara,PFMatrixList_Maize)
PFMatrixList$MA0064.1
PFMatrixList$MA1375.1

## WUS1 load
WUSmatrix1 <- rbind(A=c(11,6,78,34,7,0,70,52),
             C=c(0,79,14,1,0,65,18,13),
             G=c(1,1,1,0,0,1,0,3),
             T=c(86,12,5,64,91,32,10,30))
WUS1_1 <- PFMatrix(ID="Unknown", name="ZmWUS1_1", matrixClass="Unknown", 
                 strand="+",bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                 tags=list(), profileMatrix=WUSmatrix1)
#WUS1@name
WUSmatrix2 <- rbind(A=c(8,4,91,83,23,45,60,20,25,69),
                    C=c(1,0,0,0,2,1,1,0,30,0),
                    G=c(3,83,6,6,0,28,29,4,3,11),
                    T=c(86,11,2,9,74,24,8,75,40,18))
WUS1_2 <- PFMatrix(ID="Unknown", name="ZmWUS1_2", matrixClass="Unknown", 
                   strand="+",bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                   tags=list(), profileMatrix=WUSmatrix2)
PFMatrixList_Maize
PFMatrixList_Maize[['ZmWUS1_1']] <- WUS1_1 ## 42th
PFMatrixList_Maize[['ZmWUS1_2']] <- WUS1_2 ## 42th

#PFMatrixList_Maize[['MA0064.1']] <- PFMatrixList$MA0064.1 already Zmaize
PFMatrixList_Maize[['MA1375.1']] <- PFMatrixList$MA1375.1


PFMatrixList_All <- c(PFMatrixList_Maize,PFMatrixList_Ara)
#MOTIF HB67_WUS1_B73v5_Q30_default_finalBl_2_m0_p3_k9_c2587
#letter-probability matrix: alength= 4 w= 8 nsites= 2587
#1:A             2:C             3:G             4:T
#0.116898        0.003197        0.019711        0.860194
#0.063398        0.791174        0.016379        0.129049
#0.788296        0.143665        0.015732        0.052307
#0.343524        0.010022        0.002072        0.644382
#0.077484        0.001353        0.007036        0.914128
#0.002286        0.654206        0.013786        0.329723
#0.705835        0.186009        0.001108        0.107047
#0.523630        0.133437        0.034875        0.308058

#DE HB67_WUS1_B73v5_Q30_default_finalBl_2_m1 8.99 5 k10_c1911
#1 159(A) 33(C) 62(G) 1657(T) T
#2 86 15 1595 214 G
#3 1742 10 115 45 A
#4 1602 6 125 178 A
#5 445 40 9 1417 T
#6 877 38 536 459 A
#7 1147 36 568 160 A
#8 388 3 83 1438 T
#9 480 584 71 777 T
#10 1325 14 215 357 A
WUS1_2Temp<-data.frame(Pos1 = c(159,33,62,1657),
           Pos2=c(86,15,1595,214),
           Pos3=c(1742,10,115, 45),
           Pos4=c(1602, 6, 125, 178),
           Pos5=c(445, 40, 9, 1417),
           Pos6=c(877, 38, 536, 459),
           Pos7=c(1147, 36, 568, 160),
           Pos8=c(388, 3, 83, 1438),
           Pos9=c(480, 584, 71, 777),
           Pos10=c(1325, 14, 215, 357))
t(t(WUS1_2Temp)/colSums(WUS1_2Temp))*100

## Make BSgenome --> very annoying.... ## The chromosome should be perfectly match with input
#my_file <- read.dcf("seed_file_maizeV5.txt", fields = NULL, all = FALSE, keep.white = NULL) 
#my_file
#write.dcf(my_file , file = "seed.dcf", append = FALSE, useBytes = FALSE, indent = 0.1 * getOption("width"), width = 0.9 * getOption("width"), keep.white = NULL)
#unlink(c("BSgenome.maizeV5"), recursive = TRUE, force = TRUE)
#forgeBSgenomeDataPkg("seed.dcf") ## --> and the go to page 9 https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf


############### All motifs ############
## Find motifs
library(BSgenome.maizeV5)
rowData(ChromVarInput_counts)

ChromVarInput_counts_addGC <- addGCBias(ChromVarInput_counts, 
                            genome = BSgenome.maizeV5)
head(rowData(ChromVarInput_counts_addGC))


############ All motifs
#motif_ix <- matchMotifs(PFMatrixList_WUS, ChromVarInput_counts_addGC, 
#                        genome = BSgenome.maizeV5)

motif_ix <- matchMotifs(PFMatrixList_Maize, ChromVarInput_counts_addGC, 
                        genome = BSgenome.maizeV5)

#motif_ix <- matchMotifs(PFMatrixList_All, ChromVarInput_counts_addGC, 
#                        genome = BSgenome.maizeV5)

str(motif_ix)
str(ChromVarInput_counts)
str(motif_ix@assays@data)
dim(motif_ix@assays@data@listData$motifMatches)
Temp <- motif_ix@assays@data@listData$motifMatches
str(Temp)
head(Temp)
tail(Temp)
dim(Temp)
#Temp[,42]
#Temp@i
str(ChromVarInput_counts)
granges(ChromVarInput_counts)[Temp[,42]][500:2300]

#addGCBias() step, getBackgroundPeaks() & computeDeviations() wouldn't work without it.

dev <- computeDeviations(object = ChromVarInput_counts_addGC, annotations = motif_ix) ## It takes long when it has lots of motifs

#Error in reducer$value.cache[[as.character(idx)]] <- values : 
#  wrong args for environment subassignment
#In addition: Warning message:
#  In parallel::mccollect(wait = FALSE, timeout = 1) :
#  1 parallel job did not deliver a result

saveRDS(dev, file=paste0(SampleName,"_ChromVarDev_Maize.rds"))
#dev <- readRDS(paste0(SampleName,"_ChromVarDev.rds"))

#str(dev@elementMetadata@listData$fractionMatches)

dim(deviations(dev))
head(deviations(dev))[,c(1:10)]
rownames(deviations(dev))
str(assays(ChromVarInput_counts_addGC)$counts)
#Matrix::Matrix(assays(ChromVarInput_counts_addGC)$counts)


#### Plot Data with for loop All motifs  #########--------------------------------------------------------------------
library("here")
library(devtools)
#library(Seurat)
load_all('/home/sb14489/Socrates')

Dev <- deviations(dev)
Dev_ZScore <- deviationScores(dev)
dim(Dev_ZScore)
dim(Dev)
head(Dev)[,c(1:10)]
#######################
### Smooth Data!!!!! --> Probably to smooth data, we need tfidf--------------------------------------------------------------------
## 1) Make Markov matrix from UMAP already have

k=25
step=2
npcs=19

#head(meta_data)
#pcs <- meta_data[,c("umap1","umap2")]
head(pcs)
head(Dev)[,c(1:10)]
pcs <- pcs[rownames(meta_data),]
#pcs <- pcs[colnames(Dev),]
pcs <- pcs[,c(1:npcs)]
dim(pcs)
head(pcs)

knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx ## Calcuate the distance between cells
dim(knn.graph)
#?nn2
j <- as.numeric(x = t(x = knn.graph))
length(j)
i <- ((1:length(x = j)) - 1) %/% k + 1
edgeList = data.frame(i, j, 1)
dim(edgeList)
A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
head(A)[,c(1:10)]
dim(A)
# Smooth graph
##some questions
##1. what's meanning of A
##2. A/Matrix::rowSums(A)
##3. meanning of A%*%A
message("   * smoothing graph ...")
A = A + t(A) ##A now is symmetric
A = A / Matrix::rowSums(A)
step.size = step
if(step.size > 1){
  for(i in 1:step.size){
    message("     ~ step ",i) ## It take long more than one hour in Step3
    A = A %*% A
  }
}
saveRDS(A,file=paste0(SampleName,"_MarkovMatrix_Step2.rds"))
SampleName <- "A619"
A <- readRDS(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/NotNearPromoter/",SampleName,"_MarkovMatrix.rds"))
A <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_Markov/Bif3.MarkovMatrix.rds")

# smooth data
message("   * smoothing activity ...")
Dev_imputed <- t(A %*% t(Dev))


colnames(Dev_imputed) <- colnames(Dev)
rownames(Dev_imputed) <- rownames(Dev)
print(head(Dev_imputed[,1:10]))
dim(Dev_imputed)
str(Dev_imputed)
rownames(Dev_imputed)
colnames(Dev)[1:10]



## Ex plot 
SampleName <-"A619_Step2"

ExMotif <-Dev_imputed["ZmWUS1_1",]
meta_data$MotifDev <- c(ExMotif[rownames(meta_data)])
head(meta_data)
getwd()

ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
  geom_point(size=0.02) +
  scale_color_gradient2(low = "blue",
                        mid = "#b0d9f5",
                        high = "red")+
  theme_minimal() + ggtitle("ZmWUS1") 
ggsave(paste0(SampleName,"_ZmWUS1_1_Motif.pdf"), width=7, height=7)

ExMotif <-Dev_imputed["ZmWUS1_2",]
meta_data$MotifDev <- c(ExMotif[rownames(meta_data)])
head(meta_data)
getwd()

ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
  geom_point(size=0.02) +
  scale_color_gradient2(low = "blue",
                        mid = "#b0d9f5",
                        high = "red")+
  theme_minimal() + ggtitle("ZmWUS1_2") 
ggsave(paste0(SampleName,"_ZmWUS1_2_Motif.pdf"), width=7, height=7)

ExMotif <-Dev_imputed["MA1375.1",]
meta_data$MotifDev <- c(ExMotif[rownames(meta_data)])
head(meta_data)
ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
  geom_point(size=0.02) +
  scale_color_gradient2(low = "blue",
                        mid = "#b0d9f5",
                        high = "red")+
  theme_minimal() + ggtitle("MA1375.1") 
ggsave(paste0(SampleName,"_MA1375.1_Motif.pdf"), width=7, height=7)

ExMotif <-Dev_imputed["MA0064.1",]
meta_data$MotifDev <- c(ExMotif[rownames(meta_data)])
head(meta_data)
ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
  geom_point(size=0.02) +
  scale_color_gradient2(low = "blue",
                        mid = "#b0d9f5",
                        high = "red")+
  theme_minimal() + ggtitle("MA0064.1") 
ggsave(paste0(SampleName,"_MA0064.1_Motif.pdf"), width=7, height=7)
## Tfidf is needed ?..
## Seems like impute.activity is gene#*cell#

## Only Maize
Plotlist <- list()
for (Motif in rownames(Dev_imputed)){
  MotifRow <- Dev_imputed[Motif,]
  meta_data$MotifDev <- c(MotifRow[rownames(meta_data)])
  Plotlist[[Motif]] <- ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
    geom_point(size=0.02) +
    scale_color_gradient2(low = "blue",
                          mid = "#b0d9f5",
                          high = "red")+
    theme_minimal() + ggtitle(Motif) 
}

library(cowplot)

captured_final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
width_cal <- 6 * 5
length_cal <- (length(rownames(Dev_imputed))/6 * 5)

ggsave(paste0(SampleName,"_AllMotifs","_Maize_Imputed.pdf"), plot = captured_final_plot, 
       width = width_cal, height = length_cal, 
       units = c('in'), limitsize = FALSE,
       dpi = 300)

