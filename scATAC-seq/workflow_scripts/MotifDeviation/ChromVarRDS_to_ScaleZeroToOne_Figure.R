library(cowplot)

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
## WUS1 load  ## low cut off from GEM
WUSmatrix1 <- rbind(A=c(11,6,78,34,7,0,70,52),  ## I think it has some typos..
                    C=c(0,79,14,1,0,65,18,13),
                    G=c(1,1,1,0,0,1,0,3),
                    T=c(86,12,5,64,91,32,10,30))
WUS1_1 <- PFMatrix(ID="Unknown", name="ZmWUS1_1", matrixClass="Unknown",
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   tags=list(), profileMatrix=WUSmatrix1)
#WUS1@name

## HB67_WUS1_B73v5_Q30_qval5_finalBl_2.all.PFM_MEME.txt
WUSmatrix2 <- rbind(A=c(0,0,99,74,0,19,72,46),
                    C=c(0,0,0,1,0,0,0,12),
                    G=c(0,99,0,0,0,70,12,3),
                    T=c(99,0,0,1,99,6,7,37))
WUS1_2 <- PFMatrix(ID="Unknown", name="ZmWUS1_qval5_GEM1", matrixClass="Unknown",
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   tags=list(), profileMatrix=WUSmatrix2)
## "TCA" !
WUSmatrix3 <- rbind(A=c(12,26,23,0,0,92,53,8,5),
                    C=c(0,0,2,0,99,0,3,73,66),
                    G=c(0,0,56,0,0,0,0,3,6),
                    T=c(87,73,17,99,0,7,43,80,20))
WUS1_3 <- PFMatrix(ID="Unknown", name="ZmWUS1_qval5_GEM2", matrixClass="Unknown",
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   tags=list(), profileMatrix=WUSmatrix2)

WUSmatrix4 <- rbind(A=c(1,12,99,86,0,2,82,73,3,7,79,36,9),
                    C=c(0,0,0,8,0,0,12,15,0,38,3,32,3),
                    G=c(1,87,0,2,1,97,1,1,0,53,2,1,1),
                    T=c(97,0,0,0,98,0,2,9,95,0,14,29,85))
WUS1_4 <- PFMatrix(ID="Unknown", name="ZmWUS1_MemeChip", matrixClass="Unknown",
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   tags=list(), profileMatrix=WUSmatrix3)
#######
OC_Bif3HigherMotif1 <- rbind(A=c(87,0,6,99,1,11,4),
                             C=c(0,0,88,0,3,15,79),
                             G=c(12,0,2,0,0,68,0),
                             T=c(0,99,2,0,94,4,15))
OC_Bif3HigherMotif1_PFM <- PFMatrix(ID="Unknown", name="ATCATGC_Bif3HigherMotif1", matrixClass="Unknown",
                                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                    tags=list(), profileMatrix=OC_Bif3HigherMotif1)

OC_Bif3HigherMotif2 <- rbind(A=c(7,63,99,0,35,99,0,0,20),
                             C=c(47,36,0,0,14,0,0,0,25),
                             G=c(25,0,0,0,14,0,0,36,47),
                             T=c(20,0,0,99,35,0,99,63,7))
OC_Bif3HigherMotif2_PFM <- PFMatrix(ID="Unknown", name="SMATWATKS_Bif3HigherMotif2", matrixClass="Unknown",
                                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                    tags=list(), profileMatrix=OC_Bif3HigherMotif2)
OC_Bif3HigherMotif3 <- rbind(A=c(4,37,0,99,0,0,0,54),
                             C=c(22,0,99,0,0,0,62,18),
                             G=c(18,62,0,0,0,99,0,22),
                             T=c(54,0,0,0,99,0,37,4))
OC_Bif3HigherMotif3_PFM <- PFMatrix(ID="Unknown", name="TRCATGYA_Bif3HigherMotif3", matrixClass="Unknown",
                                    bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                    tags=list(), profileMatrix=OC_Bif3HigherMotif3)
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
PFMatrixList_Maize[['ZmWUS1_GEM1']] <- WUS1_1 ## 42th
PFMatrixList_Maize[['ZmWUS1_qval5_GEM1']] <- WUS1_2 ## 42th
PFMatrixList_Maize[['ZmWUS1_qval5_GEM2']] <- WUS1_3 ## 42th
PFMatrixList_Maize[['ZmWUS1_memeChIP']] <- WUS1_4 ## 42th

PFMatrixList_Maize[['ATCATGC']] <- OC_Bif3HigherMotif1_PFM
PFMatrixList_Maize[['SMATWATKS']] <- OC_Bif3HigherMotif2_PFM
PFMatrixList_Maize[['TRCATGYA']] <- OC_Bif3HigherMotif3_PFM
PFMatrixList_Maize[['TGAA_Arti']] <- TGAA_PFM
PFMatrixList_Maize[['TGAATGAA_Arti']] <- TGAATGAA_PFM

PFMatrixList <- c(PFMatrixList_Ara,PFMatrixList_Maize)

scale_to_minus1_1 <- function(x) {
  -1 + 2 * (x - min(x)) / (max(x) - min(x))
}

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
        scale_color_gradient2(low = "blue",
                              mid = "#c7b9bf",
                              high = "red")+
        theme_minimal() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
      ggtitle(paste0(PFMatrixList[[Motif]]@name,"_",Motif))
    }
    
    captured_final_plot <- plot_grid(plotlist = Plotlist, ncol = 6)
    width_cal <- 6 * 5
    length_cal <- (length(MotifVectors)/6 * 5)
    
    ggsave(paste0(SampleName,"_",DBName,"sAllMotifs","_Set",i,"_Imputed.pdf"), plot = captured_final_plot,
           width = width_cal, height = length_cal,
           units = c('in'), limitsize = FALSE,
           dpi = 300)
  }
}
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/ReScale_ColorChange")

## A619
A619_Maize <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif/A619_ChromVarDev_Maize.rds")
meta_data_A619 <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt")
MarkovFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_IncludingZmCLE7/A619_IncludingZmCLE7.MarkovMatrix.rds"
SampleName <- "A619"
DrawFigure_ForMotifs(A619_Maize,MarkovFile,"Maize",meta_data_A619,SampleName)
A619_Ara <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif/A619_ChromVarDev_Ara.rds")
DrawFigure_ForMotifs(A619_Ara,MarkovFile,"Ara",meta_data_A619,SampleName)

## Bif3 
Bif3_Maize <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif/Bif3_ChromVarDev_Maize.rds")
meta_data_Bif3 <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt")
MarkovFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_IncludingZmCLE7/Bif3_IncludingZmCLE7.MarkovMatrix.rds"
SampleName <- "Bif3"
DrawFigure_ForMotifs(Bif3_Maize,MarkovFile,"Maize",meta_data_Bif3,SampleName)
Bif3_Ara <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif/Bif3_ChromVarDev_Ara.rds")
DrawFigure_ForMotifs(Bif3_Ara,MarkovFile,"Ara",meta_data_Bif3,SampleName)

WUSfreqMatrix <- apply(WUSmatrix2, 2, function(x) x / sum(x))
print(WUSfreqMatrix)
WUS <- makePWM(as.matrix(WUSfreqMatrix))

pdf("ZmWUS1_logo.pdf", width = 50, height = 15)
seqLogo(WUS)
dev.off()

HDG1_MA1369Dot1 <- rbind(A=c(179,23,518,400,0,126,598,599,0,62,167),
                    C=c(56,218,67,0,0,0,0,0,0,4,291),
                    G=c(274,0,0,0,0,0,1,0,0,488,23),
                    T=c(90,358,14,199,599,473,0,0,599,45,118))
HDG1_MA1369Dot1Matrix <- apply(HDG1_MA1369Dot1, 2, function(x) x / sum(x))
print(HDG1_MA1369Dot1Matrix)
HDG1 <- makePWM(as.matrix(HDG1_MA1369Dot1Matrix))

pdf("HDG1_logo.pdf", width = 50, height = 15)
seqLogo(HDG1)
dev.off()
