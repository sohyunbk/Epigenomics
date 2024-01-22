library(devtools)
library("fgsea")
library("here")
library(devtools)
library(tidyverse)
library(Matrix)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library("optparse")
library(GenomicRanges)
library(topGO)
library(ggplot2)


option_list = list(
  make_option(c("--GOData"), type="character",
              help="GOData", metavar="character"),
  make_option(c("--WDir_fordACR"), type="character",
              help="WDir_fordACR", metavar="character"),
  make_option(c("--FileNameFix"), type="character",
              help="FileNameFix", metavar="character"),
  make_option(c("--OutPutDir"), type="character",
              help="OutPutDir", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Files
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,header=TRUE)
GOData <- "/scratch/sb14489/0.Reference/Maize_B73/GOTerm/B73_GO.RemoveBlank.out"
OutputDir <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/3.GO_from_dACR/A619_vs_Bif3_AnnV4"

## Function 1: GO list
MakeGOlist <-function(GOTerm){
  zm_func_genes <- read.table(GOData,
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

## Function2: Get nearest genes by cutoff of dACR.
Get_NearestGenes <-function(DEGInfo,FDRCutOff,logFC_PorN){
genes_ranges <- GRanges(seqnames=genes_data$V1, 
                        ranges=IRanges(start=genes_data$V2, end=genes_data$V3),
                        gene_id=genes_data$V4)

head(DEGInfo)
if (logFC_PorN == "Bif3Higher"){
SelectedPeak <- DEGInfo[DEGInfo$FDR < FDRCutOff & DEGInfo$logFC > 0,]}
else{SelectedPeak <- DEGInfo[DEGInfo$FDR < FDRCutOff & DEGInfo$logFC < 0,]}
head(SelectedPeak)
split_names <- strsplit(SelectedPeak$Peak, split = "_")
chromosomes <- sapply(split_names, function(x) x[1])
starts <- as.integer(sapply(split_names, function(x) x[2]))
ends <- as.integer(sapply(split_names, function(x) x[3]))
granges_object <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))

nearest_genes <- distanceToNearest(granges_object, genes_ranges)
distances <- mcols(nearest_genes)$distance
nearest_genes_id <- genes_ranges[subjectHits(nearest_genes)]$gene_id
return(nearest_genes_id)
}
## Function3: Fisher exact test

FisherExactTest <-function(GOList_BP,nearest_genes_id){
  Pvalues <-c()
  GOName <- c()
  for (i in c(1:length(GOList_BP))){
    ##########
    ## Fisher Exact Test
    
    #   DEG   NonDEG
    #GO 2   63
    #NonGO 1437  
    CheckedGO <- GOList_BP[[i]]  
    GONames <- names(GOList_BP)[i]
    ## it has repetitive number! 
    CheckedGO_unique <- unique(CheckedGO)
    
    GO_DEG <- length(intersect(CheckedGO_unique, nearest_genes_id))
    NonGO_DEG <- length(nearest_genes_id)-GO_DEG
    
    NonDEGs <- setdiff(genes_ranges$gene_id, nearest_genes_id)
    
    GO_NonDEG <- length(intersect(CheckedGO_unique, NonDEGs))
    NonGO_NonDEG <- length(genes_ranges$gene_id)-GO_NonDEG
    
    dat <- data.frame(
      "DEG" = c(GO_DEG, NonGO_DEG),
      "NonDEG" = c(GO_NonDEG, NonGO_NonDEG),
      row.names = c("GO", "NonGO"),
      stringsAsFactors = FALSE
    )
    
    Pvalues <- c(Pvalues,fisher.test(dat)$p.value)
    GOName <-c(GOName,GONames)
  }
  
  FDR <- p.adjust(Pvalues, method="BH")
  Result <- data.frame(GOName=GOName,Pvalues=Pvalues,FDR=FDR)
  return(Result)
}

## Function4: Draw Plot and Save the table 
DrawBarPlot_saveTable <-  function(Result, OutputDir,FileName,Color,FDRCutoff=NA){
  if (FDRCutoff==0.05){
  Sig <- Result[Result$FDR <0.05,]
  CutOFFName <- "FDR0.05"}
  else{
  Sig <- Result[Result$Pvalues <0.001,]  
  CutOFFName <- "Pvalue0.001"
  }
  #Sig <- Result[Result$Pvalues <0.001,]
  Sig_ordered <- Sig[rev(order(Sig$Pvalues)), ]
  Sig_ordered <- Sig_ordered %>% 
    arrange(match(Color, "Bif3Higher"), Color)
  Sig_ordered$logP <- -log10(Sig_ordered$Pvalues)
  Sig_ordered$GONameWithS <- paste0(Sig_ordered$Color,
                                    ".",Sig_ordered$GOName)
  Sig_ordered$GONameWithS <- factor(Sig_ordered$GONameWithS,
                                    levels = unique(Sig_ordered$GONameWithS))
  write.table(Sig_ordered, file=paste0(OutputDir,"/",FileName,"_",CutOFFName,".txt"),
              quote=F, row.names=F, col.names=T, sep="\t")
  ggplot(Sig_ordered, aes(y = GONameWithS, x = logP, fill = Color)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
    labs(y = "GO Term", x = "-log10(P-value)", title = FileName) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    scale_fill_manual(values = Color)
    ggsave(paste0(OutputDir,"/",FileName,"_",CutOFFName,".pdf"),
           width = 10, height = 1+nrow(Sig_ordered)/3, dpi = 300)
}

#### Start
GOList_BP <- MakeGOlist("BP")
GOList_CC <- MakeGOlist("CC")
GOList_MF <- MakeGOlist("MF")

NearestGenes_Bif3Higher <- Get_NearestGenes(DEGInfo,0.05,"Bif3Higher")
NearestGenes_A619Higher <-Get_NearestGenes(DEGInfo,0.05,"A619Higher")
## Save NearestGene Info!
head(NearestGenes_A619Higher)
GeneInfo <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.Curated.txt", fill = TRUE)
head(GeneInfo)
dim(GeneInfo)
ColNames <- GeneInfo[1,]
GeneInfo <- GeneInfo[-1, ]
colnames(GeneInfo) <- ColNames
head(GeneInfo)
head(GeneInfo$chr)
NearestGenes_Bif3Higher_Info <- GeneInfo[GeneInfo$gene_model %in% NearestGenes_Bif3Higher,]
write.table(NearestGenes_Bif3Higher_Info, file=paste0(OutputDir,"/NearestGenes_Bif3Higher_Info.txt"),
            quote=F, row.names=F, col.names=T, sep="\t")
NearestGenes_A619Higher_Info <- GeneInfo[GeneInfo$gene_model %in% NearestGenes_A619Higher,]
write.table(NearestGenes_A619Higher_Info, file=paste0(OutputDir,"/NearestGenes_A619Higher_Info.txt"),
            quote=F, row.names=F, col.names=T, sep="\t")
## 1) BP
Result_A619Higher <- FisherExactTest(GOList_BP,NearestGenes_A619Higher)
Result_A619Higher$Color <- "A619Higher"
Result_Bif3Higher <- FisherExactTest(GOList_BP,NearestGenes_Bif3Higher)
Result_Bif3Higher$Color <- "Bif3Higher"
Result <- rbind(Result_A619Higher,Result_Bif3Higher)
Color1 <- c("#802652","#5d850f")
DrawBarPlot_saveTable(Result, OutputDir,"BP_FDR0.05dACR",Color1,FDRCutoff=0.05)                
DrawBarPlot_saveTable(Result, OutputDir,"BP_FDR0.05dACR",Color1,FDRCutoff="NO")                

Color2 <- c("#802652","#5d850f")
Result_A619Higher_CC <- FisherExactTest(GOList_CC,NearestGenes_A619Higher)
Result_A619Higher_CC$Color <- "A619Higher"
Result_Bif3Higher_CC <- FisherExactTest(GOList_CC,NearestGenes_Bif3Higher)
Result_Bif3Higher_CC$Color <- "Bif3Higher"
Result_CC <- rbind(Result_A619Higher_CC,Result_Bif3Higher_CC)
DrawBarPlot_saveTable(Result_CC, OutputDir,"CC_FDR0.05dACR",Color2,FDRCutoff=0.05)                
DrawBarPlot_saveTable(Result_CC, OutputDir,"CC_FDR0.05dACR",Color2,FDRCutoff="NO")       

Color3 <- c("#802652","#5d850f")
Result_A619Higher_MF <- FisherExactTest(GOList_MF,NearestGenes_A619Higher)
Result_A619Higher_MF$Color <- "A619Higher"
Result_Bif3Higher_MF <- FisherExactTest(GOList_MF,NearestGenes_Bif3Higher)
Result_Bif3Higher_MF$Color <- "Bif3Higher"
Result_MF <- rbind(Result_A619Higher_MF,Result_Bif3Higher_MF)
DrawBarPlot_saveTable(Result_MF, OutputDir,"MF_FDR0.05dACR",Color3,FDRCutoff=0.05)                
DrawBarPlot_saveTable(Result_MF, OutputDir,"MF_FDR0.05dACR",Color3,FDRCutoff="NO") 



## Find the genes that are interested in .
which(names(GOList_BP) == "regulation_of_transcription,_DNA-templated.6355")
which(names(GOList_BP) == "regulation_of_flower_development.9909")
CheckedGO <- GOList_BP[[2847]]  
GONames <- names(GOList_BP)[2847]
CheckedGO_unique <- unique(CheckedGO)
intersect(CheckedGO_unique, NearestGenes_A619Higher)


CheckedGO <- GOList_BP[[3061]]  
GONames <- names(GOList_BP)[3061]
CheckedGO_unique <- unique(CheckedGO)
intersect(CheckedGO_unique, NearestGenes_Bif3Higher)
