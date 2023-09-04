library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
library("RColorBrewer")
library(gridExtra)
option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Name", metavar="character"),
  make_option(c("--PreFix"), type="character",
              help="PreFix", metavar="character"),
  make_option(c("--MinT"), type="character",
              help="MinT", metavar="character"),
  make_option(c("--MaxT"), type="character",
              help="MaxT", metavar="character"),
  make_option(c("--nPC"), type="character",
              help="nPC", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
Name <- opt$Name
PreOptions <- opt$PreFix
MinT<- opt$MinT
MaxT <- opt$MaxT
WD <- opt$WD
NumberOfPC <- opt$nPC

Ex <- function(){
  Name <- as.character("bif3_Re3")
  MinT<- as.character(0.007)
  MaxT <- as.character(0.005)
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/"
  PreOptions <- "Tn5Cut1000_Binsize500"
  #NumbeerOfWindow <- as.character(140000)
  #SVDorNMF <-as.character("SVD")
  NumberOfPC <- as.character(100)
}

##### Clean First ########################################################
print(Name)
setwd(paste0(WD,"/",Name))
obj <- readRDS(paste0(Name,"_",PreOptions,".rds"))
str(obj)
dim(obj$counts)

cell.counts <- Matrix::colSums(obj$counts)
#site.freq <- Matrix::rowMeans(obj$counts)

head(cell.counts)
#head(site.freq)

### 1) Check Tn5 Counts per cell or Cell counts per Tn5.
############ * Frequency of Tn5 per cells
df <- data.frame(counts = as.numeric(cell.counts))
ggplot(df, aes(x=counts)) +
  geom_histogram(binwidth = 500, # You can adjust binwidth as needed
                 fill="skyblue", color="black", alpha=0.7) +
  theme_minimal() +
  geom_vline(aes(xintercept=1000), color="red", linetype="dashed", linewidth=1)+
  labs(title="Frequency of Tn5 per cells",
       x="Counts", y="Frequency")
ggsave(filename = paste0(Name,"_",PreOptions,"_FrequencyofTn5perCells.pdf"), 
        width = 10, height = 7)
############ * Frequency of cells per window
NewFileName1 <-paste0(Name,"_",PreOptions,"_MinT",MinT,"_MaxT",MaxT)
x <- obj$counts
PreviousFeatureNumber <- nrow(x)
Cutoff1 <- ncol(x)*as.numeric(MinT)
dim(x)
x <- x[Matrix::rowSums(x)>(ncol(x)*as.numeric(MinT)),]
dim(x)
Cutoff2 <- quantile(Matrix::rowSums(x), c(1-as.numeric(MaxT)))
x <- x[Matrix::rowSums(x)<(quantile(Matrix::rowSums(x), c(1-as.numeric(MaxT)))),]
NewFeatureNumber <- nrow(x)
head(rowSums(obj$counts))

row_sums_data <- rowSums(obj$counts)
df <- data.frame(row_sums=row_sums_data)

plot_all <- ggplot(df, aes(x=row_sums)) +
  geom_histogram(fill="orange", color="black", binwidth=1, alpha=0.7) +  # Setting binwidth to 1 for discrete values
  labs(title = "Frequency of Tn5-Cell perWindow",
       x = "Frequency",
       y = "Count") +
  theme_minimal() +
  geom_vline(aes(xintercept=Cutoff1), color="red", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=Cutoff2), color="green", linetype="dashed", linewidth=1) +  
  annotate("text", x = Inf, y = Inf, label =
             paste0("Before window filter:", PreviousFeatureNumber, "\n",
                    "After window filter:", NewFeatureNumber,
                    "\n Cutoff1: ",Cutoff1,"\n",
                    "Cutoff2: ",Cutoff2), 
           hjust = "right", vjust = "top", size = 5)

# Then set the x limits to 0 to 20
PlotPart <- plot_all + coord_cartesian(xlim=c(0, 20))

pdf(paste0(Name,"_",NewFileName1,"_FrequencyofTn5perWindow.pdf"), width = 14, height = 7)  # Adjust width and height as needed
grid.arrange(plot_all, PlotPart, ncol=2)
dev.off()


### remove cells with less than 1,000 open peaks
### The distribution of average peak accessibilities doesnt show any clear (lower-tail) cut-offs,
#therefore, we will use the default thresholds (min.t=0.001, max.t=0.005). The arguments min.t and max.t set the minimum peak accessibility frequency and upper (99.5%) quantile cut-offs, respectively.
#Default: min.t=0.01, max.t=0.05

obj_AfterClean <- cleanData(obj, min.c=1000,
                            min.t=as.numeric(MinT), max.t=as.numeric(MaxT), verbose=T)


#########################
## 1) Normalization-tfidf
#########################
NumberOfWindow <- as.character(150000)
SVDorNMF <-as.character("SVD")


obj_new <- tfidf(obj_AfterClean)
#print("Done with normalization")

#########################
## 2) Reducing dimension
#########################
obj_MatrixDecomposition <- reduceDims(obj_new,method=SVDorNMF,n.pcs=as.numeric(NumberOfPC),
                                      num.var=as.numeric(NumberOfWindow),
                                      svd_slotName="MatrixDecomposition")

str(obj_MatrixDecomposition)
obj_UMAP <- projectUMAP(obj_MatrixDecomposition, verbose=T, k.near=50, m.dist=0.01,
                        svd_slotName="MatrixDecomposition",umap_slotName="UMAP")

##PlotDrawing functions
obj_Cluster_beforeD <- callClusters(obj_UMAP,
                                    res=0.3,
                                    verbose=T,
                                    svd_slotName="MatrixDecomposition",
                                    umap_slotName="UMAP",
                                    cluster_slotName="Clusters",
                                    cleanCluster=F,
                                    e.thresh=5)


pdf(paste0(NewFileName,"_BeforeRemovingDoublets.pdf"), width=10, height=10)
plotUMAP(obj_Cluster_beforeD, cluster_slotName="Clusters", cex=0.2)
dev.off()

saveRDS(obj_UMAP, file=paste(NewFileName,"_beforeRemovingDoublets.rds",sep=""))
#obj_UMAP <- readRDS(paste(NewFileName,"_beforeRemovingDoublets.rds",sep=""))
## 3) Remove Doublets
obj_UMAP$meta$lib_ID <- "1"
head(obj_UMAP$meta)
str(obj_UMAP)
str(obj_UMAP$rdMethod)
obj_detectDoublets <- detectDoublets(obj_UMAP, threads=10, nTrials=5,
                                     nSample=1000, rdMethod = SVDorNMF,
                                     svd_slotName="MatrixDecomposition")
dim(obj_detectDoublets$meta)
table(obj_detectDoublets$meta$doubletscore)

obj_filterDoublets <- filterDoublets(obj=obj_detectDoublets,umap_slotname = "UMAP",
                                     embedding = "UMAP",filterRatio=1.5,
                                     removeDoublets=T, libraryVar="lib_ID",
                                     verbose=TRUE)
str(obj_filterDoublets)
table(obj_filterDoublets$meta$doubletscore)

saveRDS(obj_filterDoublets, file=paste(NewFileName,".rds",sep=""))

obj_Cluster_AfterD <- callClusters(obj_filterDoublets,
                                   res=0.3,
                                   verbose=T,
                                   svd_slotName="MatrixDecomposition",
                                   umap_slotName="UMAP",
                                   cluster_slotName="Clusters",
                                   cleanCluster=F,
                                   e.thresh=5)

pdf(paste0(NewFileName,"_AfterRemovingDoublets.pdf"), width=10, height=10)
plotUMAP(obj_Cluster_AfterD, cluster_slotName="Clusters", cex=0.2)
dev.off()

## Doublets Histogram
UMAP_Table <- obj_detectDoublets$UMAP
UMAP_Table$DoubletScore <- obj_detectDoublets$meta$doubletscore
ggplot(UMAP_Table, aes(x =DoubletScore)) +
  geom_histogram()
ggsave(paste0(NewFileName,"_DoubletScoreHistogram.pdf")
       , width=12, height=10)

## Tn5 insertion score and doublet score in cluster
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAP_Table$DoubletScore), max(UMAP_Table$DoubletScore)))
ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=DoubletScore)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc
ggsave(paste0(NewFileName,"_DoubletScoreCluster.pdf")
       , width=12, height=10)


##
#obj_detectDoublets <- obj_UMAP
#UMAP_Table <- obj_detectDoublets$meta
UMAP_Table$NumberofTn5Insertion <- obj_detectDoublets$meta$log10nSites
head(UMAP_Table)
tail(UMAP_Table)
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(UMAP_Table$NumberofTn5Insertion),
                                      max(UMAP_Table$NumberofTn5Insertion)))

ggplot(UMAP_Table, aes(x=umap1, y=umap2,color=NumberofTn5Insertion)) +
  geom_point(size=0.5) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  #ylim(-0.6, 0.5)+
  #xlim(-0.6, 0.5)+
  theme_bw()+
  sc

ggsave(paste0(NewFileName,"_NumberofTn5insertionCluster.pdf")
       , width=12, height=10)
