library("here")
library(devtools)
library(Seurat)
library(stringr)
load_all('/home/sb14489/Socrates')
library(harmony)
library(symphony)
library("optparse")
library("optparse")
library(rlang)
library(ggplot2)
option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--OldRDS"), type="character",
              help="OldRDS"),
  make_option(c("--Re1"), type="character",
              help="Re1", metavar="character"),
  make_option(c("--Re2"), type="character",
              help="Re2", metavar="character"),
  make_option(c("--SampleS"), type="character",
              help="SampleS", metavar="character"),
  make_option(c("--PreOptions_forRe3Re4"), type="character",
              help="PreOptions_forRe3Re4", metavar="character"),
  make_option(c("--WD_forRe3Re4"), type="character",
              help="WD_forRe3Re4", metavar="character"),
  make_option(c("--Re3"), type="character",
              help="Re3", metavar="character"),
  make_option(c("--Re4"), type="character",
              help="Re4", metavar="character")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

WD <- opt$WD
OldRDS <- opt$OldRDS
Re1 <- opt$Re1
Re2 <- opt$Re2
SampleS <- opt$SampleS
PreOptions_forRe3Re4 <- opt$PreOptions_forRe3Re4
WD_forRe3Re4 <- opt$WD_forRe3Re4
Re3 <- opt$Re3
Re4 <- opt$Re4


Ex <- function(){
  ##OutputDir
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/4Replicates/Combined_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.rds"
  ## old data
  OldRDS <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/CombineAll/Combined_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.rds"
  Re1 <-"A619_Re1"
  Re2 <-"A619_Re2"
  ## For Additional Samples
  SampleS <- "A619"
  PreOptions_forRe3Re4 <- "Tn5Cut1000_Binsize500_MinT0.005_MaxT0.01_PC100"
  WD_forRe3Re4 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS3_FRiP4/"
  Re3 <- "A619_Re3"
  Re4 <- "A619_Re4"
}

setwd(WD)
## Make directory
if (file.exists(file.path(SampleS))){
  setwd(file.path(SampleS))
} else {
  dir.create(file.path(SampleS))
  setwd(file.path(SampleS))
}
getwd()

## 1) Combine the data and features!
## New data open : Two replicates should be in the same WD dir and naming rules.
obj_Re3 <- readRDS(paste0(WD_forRe3Re4,"/",Re3,"/",Re3,"_",PreOptions_forRe3Re4,".rds"))
obj_Re4 <- readRDS(paste0(WD_forRe3Re4,"/",Re4,"/",Re4,"_",PreOptions_forRe3Re4,".rds"))
str(obj_Re4)
head(obj_Re4$meta)
## Old data open
obj_All <- readRDS(OldRDS)
obj_Re1 <- list()
obj_Re1$meta <- subset(obj_All$meta, obj_All$meta$sampleID==Re1)
head(obj_Re1$meta)
dim((obj_Re1$meta))
obj_Re1$meta <- obj_Re1$meta[, names(obj_Re1$meta) %in% names(obj_Re4$meta)]
dim((obj_Re1$meta))
dim((obj_Re4$meta))
obj_Re1$counts <- obj_All$counts[,colnames(obj_All$counts) %in% rownames(obj_Re1$meta)]
obj_Re1$counts <- obj_Re1$counts[Matrix::rowSums(obj_Re1$counts)>0,]
str(obj_Re1)
obj_Re2 <- list()
obj_Re2$meta <- subset(obj_All$meta, obj_All$meta$sampleID==Re2)
obj_Re2$meta <- obj_Re2$meta[, names(obj_Re2$meta) %in% names(obj_Re4$meta)]
obj_Re2$counts <- obj_All$counts[,colnames(obj_All$counts) %in% rownames(obj_Re2$meta)]
obj_Re2$counts <- obj_Re2$counts[Matrix::rowSums(obj_Re2$counts)>0,]
str(obj_Re2)

SharedFeatures <- Reduce(intersect, list(rownames(obj_Re1$counts),rownames(obj_Re2$counts),
                                         rownames(obj_Re3$counts),rownames(obj_Re4$counts)))
length(SharedFeatures)
print(SharedFeatures[1:5])

#print(SharedFeatures)
files <- list(obj_Re1,obj_Re2,obj_Re3,obj_Re4)
names(files) <- c(Re1,Re2,Re3,Re4)
merged.obj <- mergeSocratesRDS(obj.list=files)
str(merged.obj)
merged.obj$counts <- merged.obj$counts[rownames(merged.obj$counts) %in% SharedFeatures,]
str(merged.obj)

##Pass filtering blacklist cause the Re3 and Re4 filtered Tn5 in blacklist in the very begining.
## 2) Harmony
####################### 
SVDorNMF <-as.character("SVD")
NumberOfPC <- as.character(300)
NumberOfWindow <- as.character(nrow(merged.obj$counts))

###########################
obj <- tfidf(merged.obj, doL2=T)
#out <-  paste0("Ref_",Prefix,"_RemoveMitoChloroChIP500bpCC")
out <-  paste0(SampleS,"_Re1234")
#saveRDS(obj, file=paste0(out,".tfidf.rds"))
#obj <- readRDS()
print("DonewithTfidf")
str(obj)
#head(obj_A619_merged$residuals)
dim(obj$residuals)
#obj_A619_merged <- readRDS(paste0(out,".tfidf.rds"))

# project with NMF -----------------------------------
obj <- reduceDims(obj,method=SVDorNMF,
                              n.pcs=as.numeric(NumberOfPC),
                              cor.max=0.7,
                              num.var=as.numeric(NumberOfWindow),
                              verbose=T,
                              scaleVar=T,
                              doSTD=F,
                              doL1=F,
                              doL2=T,
                              refit_residuals=F)

getwd()
head(obj$meta)
#- identifying variable features for clustering ...
#- keeping 200000 variable features for dimensionality reduction ...
#- reduce dimensions with SVD ... 
#Warning in irlba(t(M), nv = n.pcs) :
#  did not converge--results might be invalid!; try increasing work or maxit
#- removing components correlated to read depth...
#- normalizing reduced dimensions...

harmony_embeddings <- HarmonyMatrix(obj$PCA, meta_data=obj$meta,
                         vars_use="sampleID", do_pca=F,
                         #theta=c(3, 2),
                         sigma=0.1,
                         nclust=30,
                         max.iter.cluster=100,
                         max.iter.harmony=30,
                         return_object=F) ##return_object should be false if I want to get pca
dim(harmony_embeddings)
harmony_embeddings[seq_len(5), seq_len(5)]

obj[['HM_EMBS']] <- harmony_embeddings
str(obj)

obj_UMAP_WithHarmony <- projectUMAP(obj, verbose=T, k.near=50, 
                                    m.dist=0.01, 
                                    svd_slotName="HM_EMBS",
                                    umap_slotName="UMAP")

K <- "50"
RES <- "0.9"

#saveRDS(obj_UMAP_WithHarmony, file=paste0(out,".afterHarmony.rds"))
#obj_UMAP_WithHarmony <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/A619/A619_Tn5Cut1000_Binsize500_MinT0.005_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.afterHarmony.rds")
## I don't why but :callClusters: does not work in the terminal,.. --> figure out it takes really long in terminal.
obj_Cluster_WithHarmony <- callClusters(obj_UMAP_WithHarmony, 
                                        res=as.numeric(RES),
                                        verbose=T,
                                        k.near=as.numeric(K),
                                        svd_slotName='HM_EMBS',
                                        umap_slotName="UMAP",
                                        cluster_slotName="Clusters",
                                        cleanCluster=F,
                                        e.thresh=5)

str(obj_Cluster_WithHarmony)
str(obj_UMAP_WithHarmony)
out_final <- paste0(out,"_FeaturesN",NumberOfWindow,"_k",K,"_res",RES)
saveRDS(obj_UMAP_WithHarmony, file=paste0(out_final,".afterHarmony.rds"))

colorr <- c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#f7366d","#0bd43d",
            "#deadce","#adafde","#5703ff","#bcf089","#b376f5","#ed0505","#1010b5")
head(obj_Cluster_WithHarmony$Clusters)
All <- ggplot(obj_Cluster_WithHarmony$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) + 
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(title = paste0("Re1+R2 \n FeatureNumbers :",nrow(obj$counts),
                      "\n CellNumber: ",ncol(obj$counts)),
       x = "UMAP1",
       y = "UMAP2") 

ReplicatePlotFunction <- function(sample_id, point_color) {
  clusters_table <- subset(obj_Cluster_WithHarmony$Clusters, sampleID == sample_id)
  plot <- ggplot(clusters_table, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
    geom_point(size=0.02, color=point_color) + 
    theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    labs(title = paste0(sample_id, " : ", nrow(clusters_table)))
  return(plot)
}

# Use the function to create plots
Re1_plot <- ReplicatePlotFunction(Re1, "blue")
Re2_plot <- ReplicatePlotFunction(Re2, "red")
Re3_plot <- ReplicatePlotFunction(Re3, "purple")
Re4_plot <- ReplicatePlotFunction(Re4, "green")

## Add some plots to see the cell quality
InputMeata <- obj_Cluster_WithHarmony$Clusters
### * Tn5 log
library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(InputMeata$log10nSites),
                                      max(InputMeata$log10nSites)))
Q_Tn5 <- ggplot(InputMeata, aes(x=umap1, y=umap2, 
                color=log10nSites)) +
  geom_point(size=0.02) + 
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "LogTn5")+sc 
## * Doublets
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(InputMeata$doubletscore),
                                      max(InputMeata$doubletscore)))
Q_doubletscore <-ggplot(InputMeata, aes(x=umap1, y=umap2, 
                                        color=doubletscore)) +
  geom_point(size=0.02) + 
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "Doublet score")+sc 
## Tss ratio 
head(InputMeata)
InputMeata$rTSS <- InputMeata$tss/InputMeata$total
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(InputMeata$rTSS),
                                      max(InputMeata$rTSS)))
Q_rTSS <-ggplot(InputMeata, aes(x=umap1, y=umap2, 
                                        color=rTSS)) +
  geom_point(size=0.02) + 
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "TSS ratio")+sc 

library(gridExtra)
pdf(paste0(out_final, "_WithHarmony.pdf"), width=20, height=5)
row1 <- grid.arrange(All, Re1_plot, Re2_plot, ncol=3, widths=c(1.7, 1, 1))
row2 <- grid.arrange(Re3_plot, Re4_plot, Q_Tn5, Q_doubletscore, Q_rTSS,
                     ncol=5, widths=c(1, 1, 1, 1, 1))
grid.arrange(row1, row2, nrow=2)
dev.off()

saveRDS(obj_UMAP_WithHarmony, file=paste0(out_final,".AfterHarmony.rds"))
#obj_UMAP_WithHarmony <-readRDS("bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.rds")

head(obj_Cluster_WithHarmony$meta)
write.table(obj_Cluster_WithHarmony$Clusters, paste0(out_final,".AfterHarmony.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
head(obj_Cluster_WithHarmony$HM_EMBS)
write.table(obj_Cluster_WithHarmony$HM_EMBS, paste0(out,".AfterHarmony.PCA.txt"), quote=F, row.names=T, col.names=T, sep="\t")

####################################################################################
## Draw the barplot by Replicates 
head(obj_Cluster_WithHarmony$Clusters)
#> head(Plotdata)
#library  Fre     Ratio                         Celltype
#1 A619_Re1  393  9.451659 BundleSheath_VascularSchrenchyma
#2 A619_Re2  667  9.976069 BundleSheath_VascularSchrenchyma
#3 bif3_Re1   88  4.853833 BundleSheath_VascularSchrenchyma
getPlotDataForSampleID <- function(sampleID, META) {
  metaSubset <- META[which(META$sampleID == sampleID),]
  plotData <- data.frame(t(table(metaSubset$LouvainClusters)))
  plotData$Var1 <- sampleID
  plotData$Ratio <- (plotData$Freq / sum(plotData$Freq)) * 100
  plotData
}

# Create a function to add the position for labels
addPositionForLabels <- function(data) {
  library(plyr)
  ddply(data, .(Var1), transform, pos = cumsum(Ratio) - (0.5 * Ratio))
}
META <- obj_Cluster_WithHarmony$Clusters

sampleIDs <- c(Re1, Re2, Re3, Re4)
plotDataList <- lapply(sampleIDs, getPlotDataForSampleID, META=META)

PlotData <- do.call(rbind, plotDataList)
PlotData$Ratio <- round(PlotData$Ratio, 2)
PlotData <- addPositionForLabels(PlotData)

ClusterLevel <- levels(factor(obj_Cluster_WithHarmony$Clusters$LouvainClusters))

library(ggplot2)
p <- ggplot(PlotData, aes(fill=Var2, y=Ratio, x=Var1)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorr[1:length(ClusterLevel)]) +
  geom_text(aes(label = paste(round(Ratio, 1), "%")), position = position_stack(vjust = 0.5)) +
  theme_minimal()

ggsave(paste0(out_final, "_StackedBarplotByRe.pdf"), width=4, height=5)
