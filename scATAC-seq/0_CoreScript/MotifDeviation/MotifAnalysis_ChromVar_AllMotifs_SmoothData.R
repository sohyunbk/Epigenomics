#https://greenleaflab.github.io/chromVAR/articles/Introduction.html#example-data
# conda activate JASPAR_act
## Does chromVar ignore "+","-"strand info from JASPAR?
library(chromVAR)
library(cowplot)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
library(tidyverse)
library(stringr)
library(viridis)
library(irlba)
library(Matrix)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
#library(varistran)
library(parallel)
library(png)
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)
library(RANN)
library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(BSgenome)

library("optparse")
library(rlang)
library(ggplot2)

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Sparse"), type="character",
              help="Sparse", metavar="character"),
  make_option(c("--Meta"), type="character",
              help="Meta", metavar="character"),
  make_option(c("--pcs"), type="character",
              help="pcs", metavar="character"),
  make_option(c("--Markov"), type="character",
              help="Markov", metavar="character"),
  make_option(c("--SampleName"), type="character",
              help="SampleName", metavar="character"),
  make_option(c("--IGPeak"), type="character",
              help="IGPeak", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Sparse <- read.table(opt$Sparse,header=F)
meta_data <- read.table(opt$Meta, header=TRUE)
pcs <- read.table(opt$pcs)
MarkovFile <- opt$Markov
SampleName <- opt$SampleName
WD <- opt$WD
InterGenic_peak <- read.table(opt$IGPeak,header=F)
if (!dir.exists(WD)){
  dir.create(WD)
} else {
  print("Dir already exists!")
}
setwd(WD)

#if (SampleName == "A619"){
  #Sparse_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619/A619Intergenic500bp_A619Barcode_Tn5Count.sparse"
  #MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt"
  #pcs <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt")
  #MarkovFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_IncludingZmCLE7/A619_IncludingZmCLE7.MarkovMatrix.rds"
  #meta_data <- read.table(MetaFileA619, header=TRUE)
  #Sparse <- read.table(Sparse_A619,header=F)
  #InterGenic_peak <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619/A619.500bp_peaks_Intergenic.bed",header=F)
#} else if (SampleName == "Bif3") {
#  Sparse_Bif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/Bif3_toComPeak.sparse"
#  MetaFileBif3 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
#  pcs <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt")
#  MarkovFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_Markov/Bif3.MarkovMatrix.rds"
#  Sparse <- read.table(Sparse_Bif3,header=F)
#  meta_data <- read.table(MetaFileBif3, header=TRUE)
# InterGenic_peak <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/ComA619Bif3.unique500bpPeaks_Intergenic.bed",header=F)
#}

############################

colnames(Sparse) <- c("PeakLocus","cellID","accessability")
head(Sparse)
head(meta_data)
dim(Sparse)

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

cluster_name = "Ann_v4"
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
#str(Peak_Cell_Count)
Peak_Cell_Count_Input <- as.matrix(Peak_Cell_Count[,-1])
str(Peak_Cell_Count_Input)
#head(Peak_Cell_Count_Input)[,c(1:10)]
dim(Peak_Cell_Count_Input)

ChromVarInput_counts <- SummarizedExperiment(assays =
                                          list(counts = Peak_Cell_Count_Input),
                                        rowRanges = peaks)




saveRDS(Peak_Cell_Count, file=paste0(WD,"/",SampleName,"_Peak_perCell_Counts.rds"))
saveRDS(ChromVarInput_counts, file=paste0(WD,"/",SampleName,"_ChromVarInput.rds"))

### 4) Get motifs!! ## getMatrixSet from TFBSTools
#https://bioconductor.org/packages/devel/bioc/manuals/TFBSTools/man/TFBSTools.pdf #Page10
#conda activate MotifR !!
## Should use this file from GEM: HB67_WUS1_B73v5_Q30_default_finalBl_2.all.PFM_MEME.txt : TGAA motifs until m0, m1, m2
## It is corresponding to file:///Users/sohyun/Documents/2.SingleCellATAC/6.WUS_DAP-seqData/HB67_WUS1_B73v5_Q30_default_finalBl_LowestThreshold_Default/HB67_WUS1_B73v5_Q30_default_finalBl_outputs/HB67_WUS1_B73v5_Q30_default_finalBl_2.results.htm

#ChromVarInput_counts <- readRDS(paste0(SampleName,"_ChromVarInput.rds"))
#Peak_Cell_Count <- readRDS(paste0(SampleName,"_Peak_perCell_Counts.rds"))

str(ChromVarInput_counts)
## Remove!! chr6	181356755	181357256	1.10868807050461

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

## HB67_WUS1_B73v5_Q30_qval5_finalBl_2.all.PFM.txt

WUSmatrix2 <- rbind(A=c(13,26,3237,2412,0,622,2333,1502),
                    C=c(0,0,0,350,0,136,242,416),
                    G=c(0,3211,0,0,0,2278,416,100),
                    T=c(3224,0,0,476,3237,201,247,1221))
WUS1_2 <- PFMatrix(ID="Unknown", name="ZmWUS1_qval5_GEM1", matrixClass="Unknown",
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   tags=list(), profileMatrix=WUSmatrix2)

######### TGAATGAA
TGAA <- rbind(A=c(0,0,100,100),
              C=c(0,0,0,0),
              G=c(0,100,0,0),
              T=c(100,0,0,0))
TGAA_PFM <- PFMatrix(ID="Unknown", name="TGAA", matrixClass="Unknown",
                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                    tags=list(), profileMatrix=TGAA)
TGAATGAA <- rbind(A=c(0,0,100,100,0,0,100,100),
              C=c(0,0,0,0,0,0,0,0),
              G=c(0,100,0,0,0,100,0,0),
              T=c(100,0,0,0,100,0,0,0))
TGAATGAA_PFM <- PFMatrix(ID="Unknown", name="TGAATGAA", matrixClass="Unknown",
                          bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                        tags=list(), profileMatrix=TGAATGAA)



PFMatrixList_Maize
PFMatrixList_Maize[['ZmWUS1_qval5_GEM1']] <- WUS1_2 ## 42th
PFMatrixList_Maize[['TGAA_Arti']] <- TGAA_PFM
PFMatrixList_Maize[['TGAATGAA_Arti']] <- TGAATGAA_PFM
PFMatrixList <- c(PFMatrixList_Ara,PFMatrixList_Maize)


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
motif_ix_Maize <- matchMotifs(PFMatrixList_Maize, ChromVarInput_counts_addGC,
                    genome = BSgenome.maizeV5)
motif_ix_Arabidopsis <- matchMotifs(PFMatrixList_Ara, ChromVarInput_counts_addGC,
                        genome = BSgenome.maizeV5)

#str(motif_ix)
#str(ChromVarInput_counts)
#str(motif_ix@assays@data)
#dim(motif_ix@assays@data@listData$motifMatches)
#Temp <- motif_ix@assays@data@listData$motifMatches
#str(Temp)
#head(Temp)
#Temp@i
#str(ChromVarInput_counts)
#granges(ChromVarInput_counts)[Temp@i][1100:2000]

#addGCBias() step, getBackgroundPeaks() & computeDeviations() wouldn't work without it.

dev_Maize <- computeDeviations(object = ChromVarInput_counts_addGC, annotations = motif_ix_Maize) ## It takes long when it has lots of motifs
dev_Ara <- computeDeviations(object = ChromVarInput_counts_addGC, annotations = motif_ix_Arabidopsis) ## It takes long when it has lots of motifs

#Error in reducer$value.cache[[as.character(idx)]] <- values :
#  wrong args for environment subassignment
#In addition: Warning message:
#  In parallel::mccollect(wait = FALSE, timeout = 1) :
#  1 parallel job did not deliver a result
#setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif_ArtificalMotif")
saveRDS(dev_Maize, file=paste0(SampleName,"_ChromVarDev_Maize.rds"))
saveRDS(dev_Ara, file=paste0(SampleName,"_ChromVarDev_Ara.rds"))
#dev_Ara <- readRDS("A619_ChromVarDev_Ara.rds")
#dev_Maize <- readRDS("A619_ChromVarDev_Ara.rds")
#dev <- readRDS(paste0(SampleName,"_ChromVarDev.rds"))

#str(dev@elementMetadata@listData$fractionMatches)

#dim(deviations(dev))
#head(deviations(dev))[,c(1:10)]
#rownames(deviations(dev))
#str(assays(ChromVarInput_counts_addGC)$counts)
#Matrix::Matrix(assays(ChromVarInput_counts_addGC)$counts)

## Function for Plotting!
#dev <- dev_Maize
#DBName <- "Maize"
#MarkovFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_Markov/A619.MarkovMatrix.rds"
#MarkovFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/Past/NotNearPromoter/A619_MarkovMatrix.rds"

scale_to_minus1_1 <- function(x) {
  -1 + 2 * (x - min(x)) / (max(x) - min(x))
}

#blues <- c(brewer.pal(9, "Blues")[2:8])
greys <- brewer.pal(9, "Greys")[3:4]
#reds <- c(brewer.pal(9, "Reds")[2:8])
reds <- c("#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", 
          "#EF3B2C", "#EF3B2C", "#EF3B2C", "#CB181D", "#A50F15")
blues <- c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6",
           "#4292C6", "#4292C6", "#4292C6", "#2171B5", "#08519C")
# Intermediary colors for smoother transitions
intermediary_blue_grey <- c("#e1ebf2")
intermediary_red_grey <- c( "#ede6e9")

# Combining all colors
cols <- colorRampPalette(c(rev(blues),
                           intermediary_blue_grey,
                           greys,
                           rev(greys),
                           intermediary_red_grey,
                           reds))(100)

DrawFigure_ForMotifs <- function(dev,MarkovFile,DBName,meta_data,SampleName) {
  #dev <- dev_Ara
  Dev <- deviationScores(dev)
  
  #Dev_ZScore <- deviationScores(dev)
  #dim(Dev_ZScore)
  #dim(Dev)
  #head(Dev)[,c(1:10)]
  A <- readRDS(MarkovFile)
  Dev_ordered <- Dev[,colnames(A)]
  dim(Dev_ordered)
  dim(A)
  #head(A)[,c(1:10)]
  #head(Dev_ordered)[,c(1:10)]
  #dim(Dev)
  # smooth data
  message("   * smoothing activity ...")
  Dev_imputed <- t(A %*% t(Dev_ordered))
  colnames(Dev_imputed) <- colnames(Dev_ordered)
  rownames(Dev_imputed) <- rownames(Dev_ordered)
  print(tail(Dev_imputed[,1:10]))
  dim(Dev_imputed)
  rownames(Dev_imputed)
  Dev_Scaled_matrix <- t(apply(Dev_imputed, 1, scale_to_minus1_1))
  Dev_Scaled_dgeMatrix <- new("dgeMatrix", 
                              x = as.numeric(Dev_Scaled_matrix), 
                              Dim = dim(Dev_Scaled_matrix),
                              Dimnames = list(rownames(Dev_Scaled_matrix), colnames(Dev_Scaled_matrix)))
  Dev_imputed <- Dev_Scaled_dgeMatrix
  head(Dev_imputed)[,c(1:10)]
  NumberOfLoop <- (nrow(Dev_imputed)%/%50)+1
  for (i in c(1:NumberOfLoop)){
    if (i*50 > nrow(Dev_imputed)){
      EndLength =nrow(Dev_imputed)
    }else{EndLength = i*50}
    StartLength <- (i*50)-49
    MotifVectors <- rownames(Dev_imputed)[c(StartLength:EndLength)]
    Plotlist <- list()
    for (Motif in MotifVectors){
      MotifRow <- Dev_imputed[Motif,]
      meta_data$MotifDev <- c(MotifRow[rownames(meta_data)])
      Plotlist[[Motif]] <- ggplot(meta_data, aes(x=umap1, y=umap2, color=MotifDev)) +
        geom_point(size=0.02) +
        scale_color_gradientn(colors=cols) +
        theme_minimal() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
        ggtitle(paste0(PFMatrixList[[Motif]]@name,"_",Motif))
    }
    
    captured_final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
    width_cal <- 6 * 5
    length_cal <- (length(MotifVectors)/6 * 5)
    
    ggsave(paste0(WD,"/",SampleName,"_",DBName,"sAllMotifs","_Set",i,"_Imputed.pdf"), plot = captured_final_plot,
           width = width_cal, height = length_cal,
           units = c('in'), limitsize = FALSE,
           dpi = 300)
  }
}
DrawFigure_ForMotifs(dev_Maize,MarkovFile,"Maize",meta_data,SampleName)
print("DoneMaize")
DrawFigure_ForMotifs(dev_Ara,MarkovFile,"Ara",meta_data,SampleName)

