install.packages('devtools')
library(devtools)
install_github('welch-lab/liger')
library(rliger)
library(Seurat)
library(stringr)

#Step 1: Download the data
#First, read in your datasets. For this tutorial, we will use three matrices, which can all be downloaded at https://www.dropbox.com/sh/y9kjoum8u469nj1/AADik2b2-Qo3os2QSWXdIAbna?dl=0 .
#The transcriptomic measures (SNAREseq_RNA.RDS) is the SNARE-seq scRNA dataset (31,367 genes by 10,309 cells).
#For the shared epigenomic features (SNARE_seq_shared_chromatin_features.RDS), we create a gene-centric matrix, such that we sum of the number of accessibiltiy peaks that occur over the gene body and promoter regions for each gene. 
#For a detailed walkthough of how to generate such a matrix, please see our ’Integrating scRNA and scATAC data vignette ( http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html ). 
#The resulting matrix of gene-centric chromatin accessibility is 22,379 genes by 10,309 cells
#For the unshared epigenomic features, we binned the genome into bins of 100,000 bp, and summed the number of peaks occuring in each bin. We then filtered this matrix for all bins that overlapped with ENCODE Blacklist regions, ,genes, and promoters. Our filtered matrix (SNARE_seq_unshared_chromatin_features.RDS) is 10,309 cells by 7,437 bins.

setwd("/scratch/sb14489/4.scRNAseq/ExData")
rna = readRDS("SNAREseq_RNA.RDS")
tail(rna)[,c(1:5)]
shared_atac = readRDS("SNAREseq_chromatin_accessibility_shared.RDS")
head(shared_atac)[,c(1:5)]
unshared_atac = readRDS("SNARE_seq_unshared_chromatin_features.RDS")
tail(unshared_atac)[,c(30:38)]

rownames(rna)
rownames(shared_atac)
rownames(shared_atac)%in%rownames(rna)

str(unshared_atac)
class(unshared_atac)

liger <- createLiger(list(peaks = unshared_atac))
str(liger)
liger <- normalize(liger)
norm <- liger@norm.data$peaks
head(norm)[,c(1:5)]

se = CreateSeuratObject(norm)
vars_2000 <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
top2000 <- head(VariableFeatures(vars_2000),2000)
top2000_feats <-  norm[top2000,]   

liger <- selectGenes(liger)
liger@var.genes <- top2000
liger <- scaleNotCenter(liger)
unshared_feats = liger@scale.data$peaks

###############################################
## Step3: Preprocessing and normalization
###############################################
liger <- createLiger(list(rna = rna, atac = shared_atac))
str(liger)
head(rna)[,c(1:3)]
dim(rna)
head(shared_atac)[,c(1:3)]
dim(shared_atac)

liger <- normalize(liger)
liger <- selectGenes(liger, var.thresh = 0.1, datasets.use =1 , unshared = TRUE,  unshared.datasets = list(2), unshared.thresh= 0.2)
liger <- scaleNotCenter(liger)
peak_names <- rownames(unshared_feats)
liger@var.unshared.features[[2]] = peak_names
liger@scale.unshared.data[[2]] = t(unshared_feats)
liger <- optimizeALS(liger, k=30, use.unshared = TRUE, max_iters =30,thresh=1e-10)
liger <- quantile_norm(liger)
liger <- louvainCluster(liger)
liger <- runUMAP(liger)
str(liger)
umap_plots <-plotByDatasetAndCluster(liger, axis.labels = c("UMAP1","UMAP2"))
str(umap_plots)
head(umap_plots$data)

umap_plots[[2]]
calcAlignment(liger)
