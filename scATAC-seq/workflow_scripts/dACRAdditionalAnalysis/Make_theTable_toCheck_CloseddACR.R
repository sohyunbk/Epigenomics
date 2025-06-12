## 202208 Got it from Pablo
## 202209 Edited by Sohyun
## Dotplot for all / part of markers

# load libraries
library(dplyr)
library("edgeR")
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(tidyverse)
library(tidyr)
library(pheatmap) 
library(RColorBrewer)
library("optparse")
library(preprocessCore)
library(devtools)
library("fgsea")
library("here")
library(devtools)
library(tidyverse)
library(Matrix)
library(magrittr) # needs to be run every time you start R and want to use %>%
library("optparse")
library(GenomicRanges)
library(ggplot2)
library(edgeR)
library(preprocessCore)
library(GenomicRanges)
library(gplots)

#### 1) Find the two closest genes!
DEGFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion_NearestGENEINFO.txt"
genes_data <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed")
DEGInfo <- read.table(DEGFile,fill=TRUE,header=TRUE)
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC")
FimoWD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05A619Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/"

Fimo_TAAT <-"/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/fimo_out_2/fimo.tsv"

Fimo1 <- paste0(FimoWD,"/fimo_out_1/fimo.tsv")
Fimo3 <- paste0(FimoWD,"/fimo_out_3/fimo.tsv")
Fimo4 <- paste0(FimoWD,"/fimo_out_4/fimo.tsv")
Fimo5 <- paste0(FimoWD,"/fimo_out_5/fimo.tsv")
Fimo7 <- paste0(FimoWD,"/fimo_out_7/fimo.tsv")



Make_heatmap <- function(FimoFile,Motif,DEGInfo){

FimoOut <- read.table(FimoFile,header=TRUE)
## Find overlap.
DEGInfo$chr <- sapply(strsplit(as.character(DEGInfo$Peak), "_"), `[`, 1)
DEGInfo$start <- as.integer(sapply(strsplit(as.character(DEGInfo$Peak), "_"), `[`, 2))
DEGInfo$end <- as.integer(sapply(strsplit(as.character(DEGInfo$Peak), "_"), `[`, 3))
DEG_ranges <- GRanges(
  seqnames = DEGInfo$chr,
  ranges = IRanges(start = DEGInfo$start, end = DEGInfo$end)
)

FimoOut_ranges <- GRanges(
  seqnames = FimoOut$sequence_name,
  ranges = IRanges(start = FimoOut$start, end = FimoOut$stop)
)
overlaps <- findOverlaps(DEG_ranges, FimoOut_ranges)
overlapping_DEGInfo <- DEGInfo[queryHits(overlaps), ]
overlapping_FimoOut <- FimoOut[subjectHits(overlaps), ]
head(overlapping_DEGInfo$Peak)
overlapping_FimoOut_locus <- paste(overlapping_FimoOut$sequence_name,
                                   overlapping_FimoOut$start,
                                   overlapping_FimoOut$stop,sep="_")
OverlapTable <- data.frame(ACRLocus = overlapping_DEGInfo$Peak,overlapping_FimoOut_locus)
OverlapTable_Unique <- aggregate(overlapping_FimoOut_locus ~ ACRLocus, 
                                 data = OverlapTable, FUN = function(x) paste(x, collapse = ","))

#DEGInfo[[Motif]][DEGInfo$Peak %in% OverlapTable_Unique$ACRLocus] <- "Test"
OverlapTable_Unique_renamed <- OverlapTable_Unique
names(OverlapTable_Unique_renamed)[1] <- "Peak"
head(OverlapTable_Unique_renamed)
DEGInfo <- merge(DEGInfo, OverlapTable_Unique_renamed, by = "Peak", all.x = TRUE)
head(DEGInfo)
names(DEGInfo)[which(names(DEGInfo) == "overlapping_FimoOut_locus")] <- Motif

return(DEGInfo)
}

DEGInfo <- Make_heatmap(Fimo_TAAT,"CAATAATGC",DEGInfo)
head(DEGInfo)
DEGInfo <- Make_heatmap(Fimo1,"Fimo1_GCACAGCAGCR",DEGInfo)
DEGInfo <- Make_heatmap(Fimo3,"Fimo3_GCAGCATGC",DEGInfo)
DEGInfo <- Make_heatmap(Fimo4,"Fimo4_CGCGCCGCGCC",DEGInfo)
DEGInfo <- Make_heatmap(Fimo5,"Fimo5_YAGAGAGAGA",DEGInfo)
DEGInfo <- Make_heatmap(Fimo7,"Fimo7_GCTAGCTAGC",DEGInfo)

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/4.dfACR_TargetGene_MotifCombined")
write.table(DEGInfo,paste0("dfACR_TargetGene_MotifCombined.txt"), 
            quote=F, row.names=T, col.names=T, sep="\t")


###### Make bar plot for motif ###########
library(ggplot2)

Total_Bif3Higher <- 1437
Total_A619Higher <- 1757
  
CAATAATGC_Sum <- sum(!is.na(DEGInfo$CAATAATGC))
GCACAGCAGCR_Sum <- sum(!is.na(DEGInfo$Fimo1_GCACAGCAGCR))
GCAGCATGC_Sum <- sum(!is.na(DEGInfo$Fimo3_GCAGCATGC))
CGCGCCGCGCC_Sum <- sum(!is.na(DEGInfo$Fimo4_CGCGCCGCGCC))
YAGAGAGAGA_Sum <- sum(!is.na(DEGInfo$Fimo5_YAGAGAGAGA))
GCTAGCTAGC_Sum <- sum(!is.na(DEGInfo$Fimo7_GCTAGCTAGC))
Bif3Higher <- data.frame(Motif=c("CAATAATGC"),
                         Ratio=c((CAATAATGC_Sum/Total_Bif3Higher)*100))

Bif3Higher$Fill <- "Bif3Higher"

######


A619Higher <- data.frame(Motif=c("GCACAGCAGCR",
                                 "GCAGCATGC",
                                 "YAGAGAGAGA",
                                 "CGCGCCGCGCC",
                                 "GCTAGCTAGC"
),
Ratio=c((GCACAGCAGCR_Sum/Total_A619Higher)*100,
        (GCAGCATGC_Sum/Total_A619Higher)*100,
        (YAGAGAGAGA_Sum/Total_A619Higher)*100,
        (CGCGCCGCGCC_Sum/Total_A619Higher)*100,
        (GCTAGCTAGC_Sum/Total_A619Higher)*100
        ))

A619Higher$Fill <- "A619Higher"
A619Higher_ordered <- A619Higher[order(A619Higher$Ratio), ]
All <- rbind(A619Higher_ordered,Bif3Higher)
Color1 <- c("#4559a1","#c9ac04")

All$Motif <- factor(All$Motif,levels=c(All$Motif))
ggplot(All, aes(y = Motif, x = Ratio,fill=Fill)) +
  geom_bar(stat = "identity") +
  labs(x = "Ratio", y = "Motif") +
  ggtitle("") +
  theme_minimal()+
  xlab("% of ACRs with the motif")+
  coord_cartesian(xlim = c(0, 30))+
  scale_fill_manual(values = Color1)

ggsave("A619HigherdACR_Motif_bar_plot.pdf", width = 6, height = 4)
