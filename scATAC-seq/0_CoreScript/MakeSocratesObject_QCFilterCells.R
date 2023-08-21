library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)

## Editing the previous scripts. The thing is it should be matched with the samples using previous pipelines..

option_list = list(
  make_option(c("--WD"), type="character",
              help="WD", metavar="character"),
  make_option(c("--BinSize"), type="character",
              help="Binsize it should be 100, 500, 1000bp or peak", metavar="character"),
  make_option(c("--bed"), type="character",
              help="bed file from Tn5 insertion", metavar="character"),
  make_option(c("--ann"), type="character",
              help="GTF", metavar="character"),
  make_option(c("--chr"), type="character",
              help=".fai file", metavar="character"),
  make_option(c("--Name"), type="character",
              help="Sample file", metavar="character"),
  make_option(c("--MinTn5"), type="character",
              help="Sample file", metavar="character"),
  make_option(c("--TSS"), type="character",
            help="TSS_sd", metavar="character"),
  make_option(c("--Org"), type="character",
              help="Org ratio cutoff", metavar="character"),
    make_option(c("--FRiP"), type="character",
            help="FRiP_sd", metavar="character"),
            make_option(c("--Step"), type="character",
                    help="Step can be 'OnlyQC'", metavar="character")

  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
WD <- opt$WD
BinSize <- opt$BinSize
bed <- opt$bed
ann <- opt$ann
chr <- opt$chr
Name <- opt$Name
minimumtn5counts <- opt$MinTn5
nTSS <- as.numeric(opt$TSS)
nFRiP <- as.numeric(opt$FRiP)
nOrg <- as.numeric(opt$Org)

Example <- function(){
  Name <- as.character("bif3_Re4")
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping//"
  BinSize <- as.character("500")
  bed <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/bif3_Re4_Unique.bed"
  ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf"
  chr <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai"
  minimumtn5counts <- "1000"
}

#chr4:9929296..9936240 (6.95 Kb)
### Make wd and setwd
if (!dir.exists(WD)){
  dir.create(WD)
} else {
  print("Dir already exists!")
}

WDir <- paste0(WD,"/",Name)
if (!dir.exists(WDir)){
  dir.create(WDir)
} else {
  print("Dir already exists!")
}
setwd(WDir)

Load_Data <- function(){
obj <- loadBEDandGenomeData(bed, ann, chr)  ## Takes long
str(obj)

### Remove BlackList here : In the new sample, old sample I removed blacklist after QC.
########################=====================
## 1) Blacklist removal
########################=====================
blacklist_r <- read.table("/scratch/sb14489/3.scATAC/0.Data/BlackList/Zm.final_blaclist.final.withOutCellCycle.txt")
head(blacklist_r)
blacklist.gr <- GRanges(seqnames=as.character(blacklist_r$V1),
                        ranges=IRanges(start=as.numeric(blacklist_r$V2),
                                       end=as.numeric(blacklist_r$V3)),
                        names=as.character(blacklist_r$V4))

head(blacklist.gr)
chr.seq.lengths_load <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai")
chr.seq.lengths <- as.numeric(chr.seq.lengths_load$V2)
names(chr.seq.lengths) <- chr.seq.lengths_load$V1
intervals <- tileGenome(seqlengths=chr.seq.lengths, tilewidth=500, cut.last.tile.in.chrom=TRUE)
head(intervals)
intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),]
str(intervals)
regions <- as.data.frame(intervals)
head(regions)
regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
head(regions)

Temp <- obj$counts[rownames(obj$counts) %in% regions,]
dim(Temp)
dim(obj$counts)
obj$counts <- Temp



obj <- countRemoveOrganelle(obj,org_scaffolds=c("Mt","Pt"),remove_reads=TRUE)
## This function only add meta$MtPt and remove reads from bed file in Mtpt
## Now I even do not think that it removes the reads considering callACR function call the peaks in Pt/Mt
tail(obj$PtMt)

obj <- callACRs(obj, genomesize=2.5e9,
                shift= -50,
                extsize=100,
                fdr=0.05,
                output="bulk_peaks",
                tempdir="./macs2_temp",
                verbose=T) ## It should be run to get the Ratio of Tn5s in ACR
##Q how the peaks can be generated in MtPt although we deleted all the reads

obj <- buildMetaData(obj, tss.window=2000, verbose=TRUE,organelle_scaffolds=c("Mt","Pt"))
str(obj)

## Check the cells MtPt ratio with plot
obj$meta$pPtMt <- obj$meta$ptmt / (obj$meta$total+obj$meta$ptmt)
PlotData <- obj$meta
str(PlotData)

HighDepthnumber <- length(which(PlotData$total>=1000&PlotData$pPtMt>0.2))
LowDepthnumber <- length(which(PlotData$total<1000&PlotData$pPtMt>0.2))
PlotData$Tn5Number <- paste0("<1000 # Tn5\n (#Cells < 0.2 of rMtPt: ",
                             as.character(LowDepthnumber),")")
PlotData[which(PlotData$total>1000),]$Tn5Number <- paste0(">1000 # Tn5 \n( #Cells < 0.2 of rMtPt: ",
                                                          as.character(HighDepthnumber),")")



#str(obj)
ggplot(PlotData, aes(y=pPtMt,x=Tn5Number)) +
  geom_violin()+
  ggtitle(paste0("Total barcodeNumber :" ,as.character(length(rownames(PlotData)))))

ggsave(paste0(Name,"_OrgRatio_VioletPlot.pdf"), width=8, height=7)


##########################
##Filter Cells with high MtPt Ratio
#Cells_lowMtPt <- rownames(obj$meta)[which(obj$meta$pPtMt<0.2)]
#str(obj)
#head(obj$acr)
#obj$meta <-obj$meta[Cells_lowMtPt,]
#obj$bed <- obj$bed[obj$bed$V4%in%Cells_lowMtPt,]
#dim(obj$meta)


if (BinSize == "peak"){
  obj <- generateMatrix(obj, filtered=F,peaks=T, verbose=T,
                        organelle_scaffolds =c("Mt","Pt"))
}else{
  obj <- generateMatrix(obj, filtered=F, windows=as.numeric(BinSize),
                        peaks=F, verbose=T,organelle_scaffolds =c("Mt","Pt"))
}


## Save reds up to load data
NewFileName <- paste0(Name,"_loadData.rds")
saveRDS(obj, file=NewFileName)

}

############*************************************************************************************
############*************************************************************************************
## I removed "is cell" function
#soc.obj <- convertSparseData(obj, verbose=T)
#isCells <- dget("/home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/5_CellClustering/isCells.R")

QC <- function(){

obj <- readRDS(paste0(Name,"_loadData.rds"))

NewFileName <- paste(Name,"_Tn5Cut",minimumtn5counts,sep="")
pdf(file=paste(NewFileName,".pdf",sep=""),width=14,height=4)
obj$meta$acr <- obj$meta$acrs ##Fix the error from the variable
obj <- findCells(obj,
                 doplot=T,
                 min.cells = 100,
                 max.cells=16000,
                 set.tn5.cutoff=as.numeric(minimumtn5counts), #Override spline fitting to set minimum tn5 cout per cell.
                 min.tn5=1000, #Default:1000
                 tss.z.thresh = 2,
                 tss.min.freq = nTSS,
                 frip.min.freq = nFRiP,
                 frip.z.thresh = 3,
                 filt.org=FALSE,
                 filt.tss=TRUE,
                 filt.frip=TRUE)

dev.off()
str(obj)

obj <- convertSparseData(obj, verbose=T)

### Filtered Cells by MtPt
#obj <- readRDS(paste0(Name,"_Tn5Cut",minimumtn5counts,"_Binsize",BinSize,".rds"))
Cutoffcell<- sum(obj$meta$pPtMt < nOrg)

# Bin size control + color palette
# Create the ggplot with the desired aesthetics
p <- ggplot(obj$meta, aes(x = log10nSites, y = pPtMt)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  xlab("Tn5 integration sites per barcode (log10)") +
  ylab("Organelle Ratio") +
  theme_bw() +
  theme(legend.text = element_blank())

# Add the line across y = 0.25
p <- p + geom_hline(yintercept = nOrg, linetype = "dashed", color = "red")

# Calculate the position to place the custom text on the right top corner
x_pos <- max(obj$meta$log10nSites) - 0.1
y_pos <- max(obj$meta$pPtMt) - 0.03

# Add the custom text on the right top corner
p <- p + annotate(
  "text",
  x = x_pos, y = y_pos,
  label = paste0("# Cell = ", Cutoffcell),
  hjust = 1, vjust = 1,
  size = 4, color = "black"
)

# Save the plot as a PDF
ggsave(
  paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/",
         Name,"/",Name,"Org.pdf"),
  p,
  width = 7, height = 5,
  units = "in", dpi = 300
)
#str(obj)

obj$meta <- obj$meta[obj$meta$pPtMt < nOrg,]
str(obj)
head(obj$meta$cellID)
for( i in c(1:10)){
  print(obj$meta$cellID[i])
}
head(obj$counts[,c(1:10)])
obj$counts <- obj$counts[,colnames(obj$counts)%in%obj$meta$cellID]

NewFileName <- paste0(Name,"_Tn5Cut",minimumtn5counts,"_Binsize",BinSize,".rds")
saveRDS(obj, file=NewFileName)
writeLines(obj$meta$cellID, paste0(Name,"_FilteredCellBarcode.txt"))

}

if (opt$Step == "OnlyQC"){
  QC()
}else{
Load_Data()
QC()
}


##################
library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
library("RColorBrewer")

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
  Name <- as.character("1_A619")
  MinT<- as.character(0.01)
  MaxT <- as.character(0.05)
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/"
  PreOptions <- "Tn5Cut1000_Binsize500_Mt0.2"
  #NumbeerOfWindow <- as.character(140000)
  #SVDorNMF <-as.character("SVD")
  NumberOfPC <- as.character(100)
}

##### Clean First ########################################################
setwd(paste0(WD,"/",Name))
obj <- readRDS(paste0(Name,"_",PreOptions,".rds"))
str(obj)
dim(obj$counts)

cell.counts <- Matrix::colSums(obj$counts)
site.freq <- Matrix::rowMeans(obj$counts)

head(cell.counts)
head(site.freq)

tiff(file=paste(Name,"_",PreOptions,"_Distribution.tiff",sep=""),type="cairo")
layout(matrix(c(1:2), ncol=2))
par(mar=c(3,3,1,1))
plot(density(cell.counts), main="log10 cell counts", log="x")
abline(v=1000, col="red")
plot(density(site.freq), main="peak accessibility frequency", log="x")
dev.off()

### remove cells with less than 1,000 open peaks
### The distribution of average peak accessibilities doesnt show any clear (lower-tail) cut-offs,
#therefore, we will use the default thresholds (min.t=0.001, max.t=0.005). The arguments min.t and max.t set the minimum peak accessibility frequency and upper (99.5%) quantile cut-offs, respectively.
#Default: min.t=0.01, max.t=0.05

obj_AfterClean <- cleanData(obj, min.c=100,
                            min.t=as.numeric(MinT), max.t=as.numeric(MaxT), verbose=T)


## 1) Normalization-tfidf

NumberOfWindow <- as.character(140000)
SVDorNMF <-as.character("SVD")


obj_new <- tfidf(obj_AfterClean)
#print("Done with normalization")
## 2) Reducing dimension
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

NewFileName<-paste0(Name,"_",PreOptions,"_MinT",MinT,"_MaxT",MaxT,"_PC",NumberOfPC)

pdf(paste0(NewFileName,"_BeforeRemovingDoublets.pdf"), width=10, height=10)
plotUMAP(obj_Cluster_beforeD, cluster_slotName="Clusters", cex=0.2)
dev.off()

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
