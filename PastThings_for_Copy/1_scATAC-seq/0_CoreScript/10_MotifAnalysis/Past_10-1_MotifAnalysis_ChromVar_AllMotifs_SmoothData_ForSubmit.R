#https://greenleaflab.github.io/chromVAR/articles/Introduction.html#example-data

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

library("optparse")
library(rlang)
library(ggplot2)

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name"),
  make_option(c("--meta"), type="character",
              help="meta", metavar="character"),
  make_option(c("--geneact"), type="character",
              help="geneact", metavar="character"),
  make_option(c("--markers"), type="character",
              help="markers", metavar="character"),
  make_option(c("--pcs"), type="character",
              help="pcs", metavar="character")
);
opt_parser = OptionParser(option_list=op


SampleName <- "A619"
SampleName <- "Bif3"

if (SampleName == "A619"){
  Sparse_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/A619_CommonpeakwithBif3.sparse"
  MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
  pcs <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt")
  meta_data <- read.table(MetaFileA619, header=TRUE)
  Sparse <- read.table(Sparse_A619,header=F)
  } else if (SampleName == "Bif3") {
    Sparse_Bif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/bif3_CommonpeakwithA619.sparse"
    MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
    pcs <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt")
    Sparse <- read.table(Sparse_Bif3,header=F)
    meta_data <- read.table(MetaFileBif3, header=TRUE)
   }

############################
str(Ex)
colnames(Sparse) <- c("PeakLocus","cellID","accessability")
head(Sparse)
head(meta_data)
## Remove!! chr6	181356755	181357256	1.10868807050461 !!!!!!!!!!!!!!
dim(Sparse)
#dim(Sparse[Sparse$PeakLocus!="chr6_181356755_181357256",])
#Sparse <- Sparse[Sparse$PeakLocus!="chr6_181356755_181357256",]

#Genic_peak <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/ComA619Bif3.unique500bpPeaks_InterGenic.bed",header=F)
#InterGenic_peak <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/ComA619Bif3.unique500bpPeaks_Genic.bed",header=F)
InterGenic_peak <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_bif3_AnnV3/ComA619Bif3.unique500bpPeaks_BLRemove_Intergenic.bed",header=F)
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

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/NotNearPromoter")
SampleName

saveRDS(Peak_Cell_Count, file=paste0(SampleName,"_Peak_perCell_Counts.rds"))
saveRDS(ChromVarInput_counts, file=paste0(SampleName,"_ChromVarInput.rds"))

### 4) Get motifs!! ## getMatrixSet from TFBSTools
#https://bioconductor.org/packages/devel/bioc/manuals/TFBSTools/man/TFBSTools.pdf #Page10
#conda activate MotifR !!
## Should use this file from GEM: HB67_WUS1_B73v5_Q30_default_finalBl_2.all.PFM_MEME.txt : TGAA motifs until m0, m1, m2
## It is corresponding to file:///Users/sohyun/Documents/2.SingleCellATAC/6.WUS_DAP-seqData/HB67_WUS1_B73v5_Q30_default_finalBl_LowestThreshold_Default/HB67_WUS1_B73v5_Q30_default_finalBl_outputs/HB67_WUS1_B73v5_Q30_default_finalBl_2.results.htm

#ChromVarInput_counts <- readRDS(paste0(SampleName,"_ChromVarInput.rds"))
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
opts <- list()
#Zea mays Taxonomy ID: 4577 , Ara  Taxonomy ID: 3702
opts[["species"]] <- 4577 #NCBI tax IDs (9606).
opts[["all_versions"]] <- TRUE
PFMatrixList_Maize <- getMatrixSet(JASPAR2022, opts)
PFMatrixList <- c(PFMatrixList_Ara,PFMatrixList_Maize)
PFMatrixList$MA1835.1
PFMatrixList$MA0110.1
## WUS1 load
WUSmatrix <- rbind(A=c(11,6,78,34,7,0,70,52),
             C=c(0,79,14,1,0,65,18,13),
             G=c(1,1,1,0,0,1,0,3),
             T=c(86,12,5,64,91,32,10,30))
WUS1 <- PFMatrix(ID="Unknown", name="ZmWUS1", matrixClass="Unknown",
                 strand="+",bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                 tags=list(), profileMatrix=WUSmatrix)
WUS1@name
PFMatrixList_Maize[['ZmWUS1']] <- WUS1

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
motif_ix <- matchMotifs(WUS1, ChromVarInput_counts_addGC,
                        genome = BSgenome.maizeV5)
#motif_ix <- matchMotifs(PFMatrixList_Maize, ChromVarInput_counts_addGC,
#                        genome = BSgenome.maizeV5)
motif_ix <- matchMotifs(PFMatrixList_All, ChromVarInput_counts_addGC,
                        genome = BSgenome.maizeV5)

str(motif_ix)
str(ChromVarInput_counts)
str(motif_ix@assays@data)
str(motif_ix@assays@data@listData)
#addGCBias() step, getBackgroundPeaks() & computeDeviations() wouldn't work without it.

dev <- computeDeviations(object = ChromVarInput_counts_addGC, annotations = motif_ix) ## It takes long when it has lots of motifs

#Error in reducer$value.cache[[as.character(idx)]] <- values :
#  wrong args for environment subassignment
#In addition: Warning message:
#  In parallel::mccollect(wait = FALSE, timeout = 1) :
#  1 parallel job did not deliver a result

saveRDS(dev, file=paste0(SampleName,"_ChromVarDev.rds"))
#dev <- readRDS(paste0(SampleName,"_ChromVarDev.rds"))

dim(deviations(dev))
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
step=3
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
saveRDS(A,file=paste0(SampleName,"_MarkovMatrix.rds"))

A <- readRDS(paste0(SampleName,"_MarkovMatrix.rds"))
A <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3Test/Bif3Test.MarkovMatrix.rds")

dim(A)
# smooth data
message("   * smoothing activity ...")
Dev_imputed <- t(A %*% t(Dev))
colnames(Dev_imputed) <- colnames(Dev)
rownames(Dev_imputed) <- rownames(Dev)
print(head(Dev_imputed[,1:10]))
dim(Dev_imputed)

## Ex plot
ExMotif <-Dev_imputed["ZmWUS1",]
meta_data$MotifDev <- c(ExMotif[rownames(meta_data)])
head(meta_data)
ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
  geom_point(size=0.02) +
  scale_color_gradient2(low = "blue",
                        mid = "#b0d9f5",
                        high = "red")+
  theme_minimal() + ggtitle("ZmWUS1")
ggsave(paste0(SampleName,"_ZmWUS1_Motif.pdf"), width=7, height=7)

ExMotif <-Dev_imputed["MA1375.1",]
meta_data$MotifDev <- c(ExMotif[rownames(meta_data)])
head(meta_data)
ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
  geom_point(size=0.02) +
  scale_color_gradient2(low = "blue",
                        mid = "#b0d9f5",
                        high = "red")+
  theme_minimal() + ggtitle("MA0020.1")
ggsave(paste0(SampleName,"_MA1375.1_Motif.pdf"), width=7, height=7)
## Tfidf is needed ?..
## Seems like impute.activity is gene#*cell#

##
MotifList<- list()
MotifList[[1]] <- rownames(Dev_imputed)[c(1:150)]
MotifList[[2]] <- rownames(Dev_imputed)[c(151:300)]
MotifList[[3]] <- rownames(Dev_imputed)[c(301:450)]
MotifList[[4]] <- rownames(Dev_imputed)[c(451:length(rownames(Dev_imputed)))]

for (i in c(1:4)){
print(i)
MotifVectors <-   MotifList[[i]]
Plotlist <- list()
for (Motif in MotifVectors){
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
length_cal <- (length(MotifVectors)/6 * 5)

ggsave(paste0(SampleName,"_AllMotifs","_",i,"_Imputed.pdf"), plot = captured_final_plot,
       width = width_cal, height = length_cal,
       units = c('in'), limitsize = FALSE,
       dpi = 300)
}

## No imputation

for (i in c(1:4)){
  print(i)
  MotifVectors <-   MotifList[[i]]
  Plotlist <- list()
  for (Motif in MotifVectors){
    MotifRow <- Dev[Motif,]
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
  length_cal <- (length(MotifVectors)/6 * 5)

  ggsave(paste0(SampleName,"_AllMotifs","_",i,"_NotImputed.pdf"), plot = captured_final_plot,
         width = width_cal, height = length_cal,
         units = c('in'), limitsize = FALSE,
         dpi = 300)
}
