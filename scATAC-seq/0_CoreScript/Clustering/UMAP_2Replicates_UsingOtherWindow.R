library("here")
library(devtools)
library(Seurat)
library(stringr)
load_all('/home/sb14489/Socrates')
library(harmony)
library(symphony)
library("optparse")

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--SampleS"), type="character",
              help="Name", metavar="character"),
  make_option(c("--PreFix_name"), type="character",
              help="PreFix_name", metavar="character"),
  make_option(c("--Re1"), type="character",
              help="Re1", metavar="character"),
  make_option(c("--Re2"), type="character",
              help="Re2", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

WD <- opt$WD
SampleS <- opt$SampleS
PreOptions<- opt$PreFix_name
Re1 <- opt$Re1
Re2 <- opt$Re2

Ex <- function(){
  SampleS <- "A619"
  PreOptions <- "Tn5Cut1000_Binsize500_MinT0.007_MaxT0.005_PC100"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS3_FRiP4/"
  Re1 <- "A619_Re3"
  Re2 <- "A619_Re4"

  SampleS <- "bif3"
  PreOptions <- "Tn5Cut1000_Binsize500_MinT0.005_MaxT0.01_PC100"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS3_FRiP4/"
  Re1 <- "bif3_Re3"
  Re2 <- "bif3_Re4"
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
## Two replicates should be in the same WD dir and naming rules.
obj_Ref <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.afterHarmony.processed.rds")
str(obj_Ref)
RefWindows <- rownames(obj_Ref$counts)
length(RefWindows)

print(PreOptions)
print(paste0(WD,"/",Re1,"/",Re1,"_",PreOptions,".rds"))
obj_Re1 <- readRDS(paste0(WD,"/",Re1,"/",Re1,"_",PreOptions,".rds"))
obj_Re2 <- readRDS(paste0(WD,"/",Re2,"/",Re2,"_",PreOptions,".rds"))
length(rownames(obj_Re1$counts))
length(rownames(obj_Re2$counts))

SharedFeatures1 <- Reduce(intersect, list(rownames(obj_Re1$counts),RefWindows))
SharedFeatures2 <- Reduce(intersect, list(rownames(obj_Re2$counts),RefWindows))
SharedFeatures <- Reduce(intersect,list(SharedFeatures1,SharedFeatures2))

length(SharedFeatures1)
length(SharedFeatures2)
length(SharedFeatures)

print(SharedFeatures[1:5])
#print(SharedFeatures)
files <- list(obj_Re1,obj_Re2)
names(files) <- c(Re1,Re2)
merged.obj <- mergeSocratesRDS(obj.list=files)
str(merged.obj)
merged.obj$counts <- merged.obj$counts[rownames(merged.obj$counts) %in% SharedFeatures,]
str(merged.obj)


## 2) Harmony
#######################
SVDorNMF <-as.character("SVD")
NumberOfPC <- as.character(100)
FeatureN <- nrow(merged.obj$counts)
NumberOfWindow <- as.character(0)

###########################
obj <- tfidf(merged.obj, doL2=T)
#out <-  paste0("Ref_",Prefix,"_RemoveMitoChloroChIP500bpCC")
out <-  paste0(SampleS,"_",PreOptions,"_SameFeatureswithRef")
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
            "#deadce","#adafde","#5703ff")
head(obj_Cluster_WithHarmony$Clusters)
All <- ggplot(obj_Cluster_WithHarmony$Clusters, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02) +
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(title = paste0("Re1+R2 \n FeatureNumbers :",nrow(obj$counts),
                      "\n CellNumber: ",ncol(obj$counts)),
       x = "UMAP1",
       y = "UMAP2")

ClustersTable_Re1 <- subset(obj_Cluster_WithHarmony$Clusters, sampleID == Re1)
Re1_plot <- ggplot(ClustersTable_Re1, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02, color="blue") +
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(title = paste0("Re1 : ",nrow(ClustersTable_Re1)))
ClustersTable_Re2 <- subset(obj_Cluster_WithHarmony$Clusters, sampleID == Re2)
Re2_plot <- ggplot(ClustersTable_Re2, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
  geom_point(size=0.02, color="red") +
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(title = paste0("Re2 : ",nrow(ClustersTable_Re2)))

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
pdf(paste0(out_final,"_WithHarmony.pdf"), width=30, height=5)
grid.arrange(All, Re1_plot, Re2_plot,
             Q_Tn5, Q_doubletscore,Q_rTSS,
             ncol=6, widths=c(1.7,1,1,1,1,1))
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
META <- obj_Cluster_WithHarmony$Clusters
head(META)
ClusterLevel <- levels(factor(obj_Cluster_WithHarmony$Clusters$LouvainClusters))
RE1Meta <- META[which(META$sampleID == Re1),]
RE2Meta <- META[which(META$sampleID == Re2),]
head(RE1Meta)

RE1Plot <- data.frame(t(table(RE1Meta$LouvainClusters)))
RE1Plot$Var1 <- Re1
RE1Plot$Ratio <- (RE1Plot$Freq/sum(RE1Plot$Freq))*100
RE2Plot <- data.frame(t(table(RE2Meta$LouvainClusters)))
RE2Plot$Var1 <- Re2
RE2Plot$Ratio <- (RE2Plot$Freq/sum(RE2Plot$Freq))*100
PlotData <- rbind(RE1Plot,RE2Plot)
PlotData$Ratio <- round(PlotData$Ratio,2)
head(PlotData)
library(plyr)
PlotData <- ddply(PlotData, .(Var1),
                  transform, pos = cumsum(Ratio) - (0.5 * Ratio))


ggplot(PlotData, aes(fill=Var2, y=Ratio, x=Var1)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colorr[1:length(ClusterLevel)]) +
  geom_text(aes(label = paste(round(Ratio,1),"%")), position = position_stack(vjust =  0.5))+
  theme_minimal()
ggsave(paste0(out_final,"_StackedBarplotByRe.pdf"), width=4, height=5)
