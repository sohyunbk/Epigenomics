#library(devtools)
#library(rliger)
#library(Seurat)
#library(stringr)
#library("here")
#load_all('/home/sb14489/Socrates')
library(tidyverse)

### RNA data load 
rna <- readRDS("/scratch/sb14489/4.scRNAseq/1.DevelopmentalCellData/4.PublishedCells/AllJustCombined_Seurat.rds")
str(rna)
Temp <- head(rna@assays$RNA@counts)[,c(1:4)]
MetaData <- read.table()
head(MetaData)
RNASparse <- rna@assays$RNA@counts
NewColName <- MetaData[colnames(RNASparse),]$Celltype
#etaData[colnames(RNASparse),]$Celltype
colnames(RNASparse) <- NewColName
head(RNASparse_t)[,c(1:3)]
dim(RNASparse)
str(RNASparse)
RNASparse_t <- as.matrix(RNASparse)
te <- Matrix::colSums(RNASparse)
head(te)
RNASparse_Sum <- t(rowsum(t(RNASparse_t), colnames(RNASparse_t)))
head(RNASparse_Sum)
dim(RNASparse_Sum)

set.seed(4)
z<-matrix(sample(1:10,20, replace=T), nrow=4)
colnames(z)<-c("a","c","b","a","b")
z<-aggregate(colnames(z), data=z, sum)
t(rowsum(t(z), colnames(z)))




rownames(MetaData) <- MetaData$newBarcode
head(rna@meta.data[rownames(MetaData),])
MetaData$counts <- rna@meta.data[rownames(MetaData),]$nCount_RNA
head(MetaData)
SumTable <- MetaData %>%
  group_by(Celltype) %>%
  summarise_at(vars(), list(MeanAcc = mean))







head(rna@assays$RNA@counts)[,c(1:5)]
#typeof(rna@assays$RNA@counts)
sparse_matrix <- Matrix(rna@assays$RNA@counts)
str(sparse_matrix)

read_delim(sparse_matrix)
#Temp <- head(sparse_matrix)[,c(1:3)]
raw_cpm_counts_all_genes <- read_delim(sparse_matrix, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
  dplyr::mutate(cellID = barcode)  %>%
  dplyr::mutate(geneID = gene_name)
head(raw_cpm_counts_all_genes)

Temp %>%
  as.tibble(rownames="gene_name") %>%
  gather()
 
Temp <- data.frame(Temp)
Temp %>% 
  gather(key = "year", value = "cases", 2:3)

m <- matrix(1:3, nrow = 3, dimnames = list(c("X","Y","Z"), c("A")))

newtib <- m %>%
  as.tibble(rownames = "group1") %>%
  gather('A', key = "group2", value = "value")

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










