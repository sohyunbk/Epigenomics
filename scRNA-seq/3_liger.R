library(devtools)
library(rliger)
library(Seurat)
library(stringr)
library("here")
load_all('/home/sb14489/Socrates')

### RNA data load 
rna <- readRDS("/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/4.PublishedCells/AllJustCombined_Seurat.rds")
head(rna@assays$RNA@counts)[,c(1:5)]
rna <- rna@assays$RNA@counts


### ATAC data load 
## Bin --> normalized OK ( but the bin should not include the bins within genes!!!)
## GA Matrix

obj <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.afterHarmony.processed.rds")
str(obj)
head(obj$residuals)[,c(1:10)]
Regions <- rownames(obj$residuals)
Temp <- head(Regions)

vSeqNames <- sapply(strsplit(Regions, "_"), "[", 1)
vStart <- sapply(strsplit(Regions, "_"), "[", 2)
vEnd <- sapply(strsplit(Regions, "_"), "[", 3)

vStart <- as.numeric(vStart)
vEnd <- as.numeric(vEnd)

Feature_Grange <-  GRanges(seqnames = vSeqNames,
                         ranges = IRanges(start = vStart,
                                          end = vEnd,
                                          names = Regions))
  
ann <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf.gtf"
gff <- suppressWarnings(suppressMessages(makeTxDbFromGFF(as.character(ann), format="gtf", dbxrefTag="Parent")))
Ann <- genes(gff)
nRange = 500
Start_new <- c((ranges(Ann)+nRange)@start)-1 ## To convert to the bed format coordinate
Width_new <- c((ranges(Ann)+nRange)@width)+1 ## To convert to the bed format coordinate
BroadRange_Ann <- GRanges(seqnames=Ann@seqnames,
                          ranges= 
                            IRanges(start=Start_new,
                                    width=Width_new,
                                    names=names(Ann)))
#obj$ann_broad <- BroadRange_Ann

## Find overlap
hits_Within <- findOverlaps(Feature_Grange,  BroadRange_Ann,minoverlap=1,
                            type=c("within"),select="all",ignore.strand = TRUE)

Intersect <- paste(names(BroadRange_Ann)[hits_Within@to],
                   names(Feature_Grange)[hits_Within@from],sep="/_Com_/")
head(Intersect)
Intersect <- table(Intersect)

Intersect<- data.frame(gene_name = as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                split="/_Com_/"), "[", 1)),
                       barcode= as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                             split="/_Com_/"), "[", 2)),
                       accessability = as.character(Intersect))
head(Intersect)

dim(obj$residuals)
head(obj$residuals)[,c(1:10)]
unshared_atac <- obj$residuals[!rownames(obj$residuals)%in%Intersect$barcode,]
dim(unshared_atac)
rownames(unshared_atac) <- gsub("_","-",rownames(unshared_atac))

GA <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/GA_A619.txt",header=TRUE)
library(tidyverse)
shared_atac <- GA %>% 
  spread(key = gene_name, value = accessability,fill=0)
dim(shared_atac)
rownames(shared_atac) <- shared_atac$barcode
shared_atac <- shared_atac[,-c(1)]
head(shared_atac)[,c(1:5)]

shared_atac_t <- t(shared_atac)
head(shared_atac_t)[,c(1:5)]
head(unshared_atac)[,c(1:5)]

shared_atac <- shared_atac_t
##############################################################
rownames(shared_atac)%in%rownames(rna)

### Selecting the unshared features
liger <- createLiger(list(peaks = unshared_atac))
str(liger)
## Dont' need to normalize
#norm <- liger@norm.data$peaks
liger@norm.data <- liger@raw.data
norm <- liger@norm.data$peaks

se = CreateSeuratObject(norm)
vars_2000 <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
top2000 <- head(VariableFeatures(vars_2000),2000)
top2000_feats <-  norm[top2000,]   
head(top2000_feats)
#?selectGenes
liger <- selectGenes(liger,var.thresh=0)
liger@var.genes <- top2000
liger <- scaleNotCenter(liger)
unshared_feats = liger@scale.data$peaks
str(unshared_feats)
dim(unshared_feats)
dim(rna)
head(rna)[,c(1:4)]
head(shared_atac)[,c(1:4)]


## Combine two
dim(rna)
length(rownames(rna))
dim(shared_atac)

liger <- createLiger(list(rna = rna, atac = shared_atac))
str(liger)
liger <- normalize(liger)










