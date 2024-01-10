library(devtools)
library("fgsea")
library("here")
library(devtools)
library(tidyverse)
library(Matrix)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library("optparse")


option_list = list(
  make_option(c("--GOData"), type="character",
              help="GOData", metavar="character"),
  make_option(c("--WDir_forDenovo"), type="character",
              help="WDir_forDenovo", metavar="character"),
  make_option(c("--FileNameFix"), type="character",
              help="FileNameFix", metavar="character"),
  make_option(c("--OutPutDir"), type="character",
              help="OutPutDir", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#"/scratch/sb14489/0.Reference/Maize_B73/GOTerm/B73_GO.RemoveBlank.out"

MakeGOlist <-function(GOTerm){
  zm_func_genes <- read.table(opt$GOData,
                              header=TRUE)
  BP <- zm_func_genes[which(zm_func_genes$ontology==GOTerm),]
  #head(BP)
  #dim(BP)
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
  return(BP_list)}

GettingGOResult <-function(BP_list,FileName,GOTerm){
  #str(BP_list)
  #CellName <- "XylemParenchyma"
  #DEDW <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV4/"
  DEDW <- opt$WDir_fordACR
  DESeq2Result <- read.table(paste0(DEDW,FileName),
                             header=TRUE)
  #head(DESeq2Result)
  DESeq2Result_filtered <- DESeq2Result[which(DESeq2Result$padj < 0.05 &DESeq2Result$log2FoldChange >0.1),]
  zm_bs_denovo.selected <- DESeq2Result_filtered  %>% 
    dplyr::select(gene_name, stat)
  #head(zm_bs_denovo.selected)
  ranks <- deframe(zm_bs_denovo.selected)
  ranks_ordered <- ranks[order(ranks)]
  #head(ranks_ordered)
  #length(ranks_ordered)
  # probalby sort ranks
  
  # stats : Named vector of gene-level stats.
  Test <- names(head(ranks_ordered))
  #https://code.bioconductor.org/browse/fgsea/RELEASE_3_18/
  fgseaRes <- fgsea(pathways=BP_list, stats=ranks_ordered, 
                    nperm=10000,
                    minSize = 10, maxSize = 600)
  #str(fgseaRes)  
  #head(fgseaRes)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(NES)
  
  Output <- fgseaResTidy[order(fgseaResTidy$padj,decreasing=FALSE,na.last=TRUE),]
  #head(Output)
  OutputTable <- data.frame(Output)
  #colnames(OutputTable)
  OutputTable <- subset(OutputTable, select = -c(leadingEdge))
  #str(OutputTable)
  #head(OutputTable)
  setwd(opt$OutPutDir)
  CellName <-  sub(FileEnd, "", FileName)
  write.table(OutputTable, file=paste0(CellName,".",GOTerm,".txt"),
              quote=F, row.names=F, col.names=T, sep="\t")
}


## Start here 

GOName <- c("BP","CC","MF")

DEDW <- opt$WDir_forDenovo
FileEnd <- opt$FileNameFix
#FileEnd <- "_deseq_2_results.tsv"
pattern_string <- paste0(FileEnd, "$")
files <- list.files(path = DEDW, pattern = pattern_string, full.names = FALSE)

for (File in files){
  GOList_BP<- MakeGOlist("BP")
  GettingGOResult(GOList_BP,File,"BP")
  GOList_CC<- MakeGOlist("CC")
  GettingGOResult(GOList_CC,File,"CC")
  GOList_MF<- MakeGOlist("MF")
  GettingGOResult(GOList_MF,File,"MF")
  }
