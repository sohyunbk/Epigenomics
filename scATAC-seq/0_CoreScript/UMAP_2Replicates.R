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
  make_option(c("--PreFix"), type="character",
              help="PreFix", metavar="character"),
  make_option(c("--Re1"), type="character",
              help="Re1", metavar="character"),
  make_option(c("--Re2"), type="character",
              help="Re2", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

WD <- opt$WD
SampleS <- opt$SampleS
PreFix<- opt$PreFix
Re1 <- opt$Re1
Re2 <- opt$Re2

Ex <- function(){
  SampleS <- "A619"
  Prefix <- "Tn5Cut1000_Binsize500_MinT0.01_MaxT0.05_PC100"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/"
  Re1 <- "A619_Re3"
  Re2 <- "A619_Re4"
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
print(Prefix)
print(paste0(WD,"/",Re1,"/",Re1,"_",Prefix,".rds"))
obj_Re1 <- readRDS(paste0(WD,"/",Re1,"/",Re1,"_",Prefix,".rds"))
obj_Re2 <- readRDS(paste0(WD,"/",Re2,"/",Re2,"_",Prefix,".rds"))
length(rownames(obj_Re1$counts))
SharedFeatures <- Reduce(intersect, list(rownames(obj_Re1$counts),rownames(obj_Re2$counts)))
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
NumberOfPC <- as.character(300)
NumbeerOfWindow <- as.character(0)

###########################
obj <- tfidf(merged.obj, doL2=T)
#out <-  paste0("Ref_",Prefix,"_RemoveMitoChloroChIP500bpCC")
out <-  paste0(SampleS,"_",Prefix)
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
                              num.var=as.numeric(NumbeerOfWindow),
                              verbose=T,
                              scaleVar=T,
                              doSTD=F,
                              doL1=F,
                              doL2=T,
                              refit_residuals=F)

getwd()
head(obj$meta)
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
out_final <- paste0(out,"_FeaturesN",NumbeerOfWindow,"_k",K,"_res",RES)

pdf(paste0(out_final,"_WithHarmony.pdf"), width=10, height=10)
plotUMAP(obj_Cluster_WithHarmony, cluster_slotName="Clusters", cex=0.2)
dev.off()

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

library(gridExtra)
pdf(paste0(out_final,"_WithHarmony.pdf"), width=15, height=5) 
grid.arrange(All, Re1_plot, Re2_plot, ncol=3, widths=c(1.7,1,1))  
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
ClusterLevel <- levels(factor(obj_Cluster_WithHarmony$Clusters$LouvainClusters))
RE1Meta <- META[which(META$library == paste0(SampleS,"_Re1")),]
RE2Meta <- META[which(META$library == paste0(SampleS,"_Re2")),]
head(RE1Meta)

RE1Plot <- data.frame(t(table(RE1Meta$LouvainClusters)))
RE1Plot$Var1 <- "Re1"
RE1Plot$Ratio <- (RE1Plot$Freq/sum(RE1Plot$Freq))*100
RE2Plot <- data.frame(t(table(RE2Meta$LouvainClusters)))
RE2Plot$Var1 <- "Re2"
RE2Plot$Ratio <- (RE2Plot$Freq/sum(RE2Plot$Freq))*100
PlotData <- rbind(RE1Plot,RE2Plot)
PlotData$Ratio <- round(PlotData$Ratio,2)

library(plyr)
PlotData <- ddply(PlotData, .(Var1),
                  transform, pos = cumsum(Ratio) - (0.5 * Ratio))


ggplot(PlotData, aes(fill=Var2, y=Ratio, x=Var1)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colorr[1:length(ClusterLevel)]) +
  geom_text(aes(label = paste(round(Ratio,1),"%")), position = position_stack(vjust =  0.5))+
  theme_minimal()
ggsave(paste0(out_final,"_StackedBarplotByRe.pdf"), width=4, height=5)


