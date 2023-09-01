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
head(obj$bed)
### Remove BlackList here : In the new sample, old sample I removed blacklist after QC.
########################=====================
## 1) Blacklist removal --> It should be after binning!!!!
########################=====================
blacklist_r <- read.table("/scratch/sb14489/3.scATAC/0.Data/BlackList/Zm.final_blaclist.final.withOutCellCycle.txt")
head(blacklist_r)
blacklist.gr <- GRanges(seqnames=as.character(blacklist_r$V1),
                        ranges=IRanges(start=as.numeric(blacklist_r$V2),
                                       end=as.numeric(blacklist_r$V3)),
                        names=as.character(blacklist_r$V4))

head(blacklist.gr)
Bed_granges <- GRanges(seqnames=as.character(obj$bed$V1),
                       ranges=IRanges(start=as.numeric(obj$bed$V2),
                                      end=as.numeric(obj$bed$V3)),
                       names=as.character(obj$bed$V4),
                       strand=as.character(obj$bed$V5))


Bed_NotBlacklist <- Bed_granges[-queryHits(findOverlaps(Bed_granges, blacklist.gr, type="any")),]
str(Bed_NotBlacklist)

regions <- as.data.frame(Bed_NotBlacklist)
dim(regions)
dim(obj$bed)

head(obj$bed)
head(regions)

regions_new <- regions %>%
  select(V1 = seqnames, V2 = start, V3 = end, V4 = names, V5 = strand)
head(regions_new)

obj$bed  <- regions_new

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
