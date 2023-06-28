library(devtools)
library("fgsea")
library("here")
library(devtools)
library(tidyverse)
library(Matrix)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)


zm_func_genes <- read.table("/scratch/sb14489/0.Reference/Maize_B73/GOTerm/B73_GO.RemoveBlank.out",header=TRUE)
#zm_func_genes <- read.table("e_downloads/Zm-B73-REFERENCE-NAM-5.0/biological_process.gmt.txt")
head(zm_func_genes)
BP <- zm_func_genes[which(zm_func_genes$ontology=="BP"),]
head(BP)
dim(BP)
BP$gene_id <- sapply(strsplit(as.character(BP$qpid),'_'), "[", 1) 
BP$Name <- paste0(BP$desc,".",BP$goid)
BP_list <- list()
BPNames <- levels(factor(BP$Name))
length(BPNames)
i <- BPNames[2]
for(i in BPNames){
  Sub <- BP[which(BP$Name == i),]
  #head(Sub)
  BP_list[[i]] <- Sub$gene_id
}
str(BP_list)

zm_functional_annotation <- convertTerms(zm_func_genes)
head(zm_functional_annotation)
str(zm_functional_annotation)

zm_bs_denovo <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV3/Ref_AnnV3.1_NA.upregulated_genes.deseq2_output.tsv",
                           header = TRUE)

TEMP <- read.table(paste0(DEDW,"Ref_AnnV3.ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma.upregulated_genes.deseq2_output.tsv"),header=TRUE)
tail(TEMP)
dim(TEMP)

DEDW = "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV3/"
DESeq2Result <- read.table(paste0(DEDW,
  "Ref_AnnV3.1_NA_deseq_2_results.tsv"),
  header=TRUE)
head(DESeq2Result)
zm_bs_denovo <- DESeq2Result[which(DESeq2Result$padj < 0.05 &DESeq2Result$log2FoldChange >0.1),]
dim(zm_bs_denovo)

head(zm_bs_denovo)
read_table_run_FGSEA <- function(de_seq2_input, location, functional_annotations){
  
  input_file_name <- here(location, de_seq2_input)
  zm_bs_denovo <- read.table(input_file_name,
                             header = TRUE)
  
  zm_bs_denovo.selected <- zm_bs_denovo  %>% 
    dplyr::select(gene_name, stat)
  #head(zm_bs_denovo.selected)
  ranks <- deframe(zm_bs_denovo.selected)
  ranks_ordered <- ranks[order(ranks)]
  head(ranks_ordered)
  length(ranks_ordered)
  # probalby sort ranks
  
  fgseaRes <- fgsea(pathways=BP_list, stats=ranks_ordered, 
                    nperm=10000,
                    minSize = 10, maxSize = 600)
  #str(fgseaRes)  
  #head(fgseaRes)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(NES)
  
  #head(fgseaRes)
  fgseaResTidy$celltype <- de_seq2_input
  Temp <- fgseaResTidy[order(fgseaResTidy$padj,decreasing=FALSE,na.last=TRUE),]
  head(Temp)
  Temp[c(1:20),]$pathway
  return(fgseaResTidy)
}

