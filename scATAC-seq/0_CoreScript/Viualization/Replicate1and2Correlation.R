### This script is to draw the correlation plot as bulk
### It's from bed to plot
library(GenomicRanges)
library(ggplot2)
library(dplyr)

## Should it be before QC or after QC? Let's do before QC.-from all bed file.

Replicate1_Summit <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619/macs2_temp/bulk_peaks_summits.bed"
Replicate2_Summit <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619_2/macs2_temp/bulk_peaks_summits.bed"
SummitFile <- Replicate2_Summit
ReadSummit <- function(SummitFile){

Summit <- read.table(SummitFile)
#head(Summit)
Summit$V6 <- (Summit$V5/sum(Summit$V5))*1000000 ## it's CPM normalization
Summit_chr <- Summit[grepl("^chr", Summit$V1), ]
Summit_Grange <-  GRanges(seqnames = Summit_chr$V1,
                       ranges = IRanges(start = Summit_chr$V2,
                                        end = Summit_chr$V3,
                                        names = Summit_chr$V6))
nRange = 250
Start_new <- c((ranges(Summit_Grange)+nRange)@start) ## it minus nRange from the start : -1 means to convert to the bed format coordinate
#head(Start_new)
Width_new <- c((ranges(Summit_Grange)+nRange)@width)-1 ## it add the end range with nRange : +1 means to convert to the bed format coordinate
Peak500bp_Grange <- GRanges(seqnames=Summit_Grange@seqnames,
                          ranges= 
                            IRanges(start=Start_new,
                                    width=Width_new,
                                    names=names(Summit_Grange)))
negative_indices <- which(start(Peak500bp_Grange) < 0 | end(Peak500bp_Grange) < 0)
Peak500bp_Grange <- Peak500bp_Grange[-negative_indices]
return(Peak500bp_Grange)
}

### 1) Combine peaks from two replicates
Peak500bp_Re1 <- ReadSummit(Replicate1_Summit)
Peak500bp_Re2 <- ReadSummit(Replicate2_Summit)
 
# 2)  Finding overlaps between Re1 and Re2
overlaps <- findOverlaps(Peak500bp_Re1, Peak500bp_Re2)

overlaps_Re1 <- Peak500bp_Re1[queryHits(overlaps)]
overlaps_Re2 <- Peak500bp_Re2[subjectHits(overlaps)]
df_Re1 <- as.data.frame(overlaps_Re1, row.names=NULL)
df_Re2 <- as.data.frame(overlaps_Re2, row.names=NULL)
df_Re1$origin <- "Re1"
df_Re2$origin <- "Re2"
df_Re1$num_value <- as.numeric(names(overlaps_Re1))
df_Re2$num_value <- as.numeric(names(overlaps_Re2))

head(df_Re1)
head(df_Re2)

df_result <- df_Re1
head(df_result)
nrow(df_result)
# Loop through each row
for (i in 1:nrow(df_Re1)) {
  # Compare num_value and select the row with the higher value
  if (df_Re2$num_value[i] > df_Re1$num_value[i]) {
    df_result[i, ] <- df_Re2[i, ]
  }
}
head(df_result)
dim(df_result)

### 3) There can be still overlap because we used summit here.
OverlapGRange <- GRanges(seqnames = df_result$seqnames,
        ranges = IRanges(start = df_result$start,
                         end = df_result$end,
                         names = df_result$num_value))

OverlapforOverlap <- findOverlaps(OverlapGRange, OverlapGRange)
#nRow <-nrow(df_result)+1
#max(queryHits(OverlapforOverlap))
FinalTable <- data.frame()
for(nRow in 1:nrow(df_result)){
OverlappedRowNumber <- subjectHits(OverlapforOverlap)[which(queryHits(OverlapforOverlap)==nRow)]
TempTable <- df_result[OverlappedRowNumber,]
highest_value_row <- TempTable[which.max(TempTable$num_value), ]
FinalTable <- rbind(FinalTable,highest_value_row)
}
