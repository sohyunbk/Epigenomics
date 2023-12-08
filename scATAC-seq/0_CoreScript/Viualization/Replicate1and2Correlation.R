### This script is to draw the correlation plot as bulk
### It's from bed to plot
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(preprocessCore)
library(viridis)
library(stringr)
library("optparse")

## Should it be before QC or after QC? Let's do before QC.-from all bed file.

option_list = list(
  make_option(c("--Re1_Summit"), type="character",
              help="Re1_Summit", metavar="character"),
  make_option(c("--Re2_Summit"), type="character",
              help="Re2_Summit"),
  make_option(c("--Re1_AllReads"), type="character",
              help="Re1_AllReads", metavar="character"),
  make_option(c("--Re2_AllReads"), type="character",
              help="Re2_AllReads"),
  make_option(c("--OutFileName"), type="character",
              help="OutFileName", metavar="character"),
  make_option(c("--OutPath"), type="character",
              help="OutPath", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#Replicate1_Summit <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619/macs2_temp/bulk_peaks_summits.bed"
#Replicate2_Summit <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619_2/macs2_temp/bulk_peaks_summits.bed"
#Replicate1_AllReads <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/Test.bed"
#Replicate2_AllReads <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_2_Unique.bed"

Replicate1_Summit <- opt$Re1_Summit
Replicate2_Summit <- opt$Re2_Summit
Replicate1_AllReads  <- opt$Re1_AllReads
Replicate2_AllReads <- opt$Re2_AllReads
OutFilePath <- opt$OutPath
OutFileName <- opt$OutFileName


## Function1: Read Summit

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


## Function 2: Getting common peak
GettingCommonPeak <- function(Peak500bp_Re1,Peak500bp_Re2){
  
overlaps <- findOverlaps(Peak500bp_Re1, Peak500bp_Re2)

overlaps_Re1 <- Peak500bp_Re1[queryHits(overlaps)]
overlaps_Re2 <- Peak500bp_Re2[subjectHits(overlaps)]
df_Re1 <- as.data.frame(overlaps_Re1, row.names=NULL)
df_Re2 <- as.data.frame(overlaps_Re2, row.names=NULL)
df_Re1$origin <- "Re1"
df_Re2$origin <- "Re2"
df_Re1$num_value <- as.numeric(names(overlaps_Re1))
df_Re2$num_value <- as.numeric(names(overlaps_Re2))

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
return(FinalTable)
}


ReadSummit <- function(Tn5File,CommonPeak_GRange_unique){
print("Read Tn5 file")
AllTn5 <- read.table(Tn5File)
print("Done:Read Tn5 file")
#chr1	51	52	CB:Z:TAGACTGGTGCAAGAC-1_A619	+
Tn5_GRange <- GRanges(seqnames = AllTn5$V1,
                      ranges = IRanges(start = AllTn5$V2,
                                       end = AllTn5$V3,
                                       names = AllTn5$V4))
# Find overlaps
print('Getting overlap between commonpeak and Tn5')
overlaps_final <- findOverlaps(Tn5_GRange, CommonPeak_GRange_unique)
print('Done: getting overlap between commonpeak and Tn5')

# Create a data frame to store the results
overlap_counts <- data.frame(seqnames = seqnames(CommonPeak_GRange_unique),
                             start = start(CommonPeak_GRange_unique),
                             end = end(CommonPeak_GRange_unique),
                             count = integer(length(CommonPeak_GRange_unique)))
# Count the overlaps for each range in CommonPeak_GRange
for (i in seq_along(CommonPeak_GRange_unique)) {
  overlap_counts$count[i] <- sum(subjectHits(overlaps_final) == i)
}
return(overlap_counts)
}

### 1) Combine peaks from two replicates
Peak500bp_Re1 <- ReadSummit(Replicate1_Summit)
Peak500bp_Re2 <- ReadSummit(Replicate2_Summit)
#### 2)  Finding overlaps between Re1 and Re2
CommonPeakOutfile <- paste0(OutFilePath,OutFileName,"_BulkCommonPeak.bed")
if (file.exists(CommonPeakOutfile)) {
  FinalTable <- read.table(CommonPeakOutfile,header=TRUE)
}else{
FinalTable <- GettingCommonPeak(Peak500bp_Re1,Peak500bp_Re2)
write.table(FinalTable,CommonPeakOutfile, quote=F, row.names=F, col.names=T, sep="\t")
}
CommonPeak_GRange <- GRanges(seqnames = FinalTable$seqnames,
                             ranges = IRanges(start = FinalTable$start,
                                              end = FinalTable$end,
                                              names = FinalTable$num_value))
print("CommonPeak:")
head(CommonPeak_GRange)
CommonPeak_GRange_unique <- unique(CommonPeak_GRange)
###### 3) getting the read depth for the common peaks
if (file.exists(paste0(OutFilePath,OutFileName,"_Tn5CountToCommonPeak_Re1.bed"))) {
  overlap_counts_Re1 <- read.table(paste0(OutFilePath,OutFileName,"_Tn5CountToCommonPeak_Re1.bed"),header=T)
}else{
  overlap_counts_Re1 <- ReadSummit(Replicate1_AllReads,CommonPeak_GRange_unique)
  write.table(overlap_counts_Re1,paste0(OutFilePath,OutFileName,"_Tn5CountToCommonPeak_Re1.bed"), quote=F, row.names=F, col.names=T, sep="\t")
}

if (file.exists(paste0(OutFilePath,OutFileName,"_Tn5CountToCommonPeak_Re2.bed"))) {
  overlap_counts_Re2 <- read.table(paste0(OutFilePath,OutFileName,"_Tn5CountToCommonPeak_Re2.bed"),header=T)
} else{
  overlap_counts_Re2 <- ReadSummit(Replicate2_AllReads,CommonPeak_GRange_unique)
  write.table(overlap_counts_Re2,paste0(OutFilePath,OutFileName,"_Tn5CountToCommonPeak_Re2.bed"), quote=F, row.names=F, col.names=T, sep="\t")
}

### 4) Draw Correlation plot

CorrTable <- data.frame(peak=paste(overlap_counts_Re1$seqnames,
                                   overlap_counts_Re1$start,
                                   overlap_counts_Re1$end,
                                   sep="_"),
                        Re1_count=overlap_counts_Re1$count,
                        Re2_count=overlap_counts_Re2$count)

rownames(CorrTable) <- CorrTable$peak
CorrTable <- CorrTable[,-1]
# CPM Normalization
quantile_normalized <- normalize.quantiles(as.matrix(CorrTable))
quantile_normalized_df <- as.data.frame(quantile_normalized)
quantile_normalized_df$logRe1 <- log10(quantile_normalized_df$V1)
quantile_normalized_df$logRe2 <- log10(quantile_normalized_df$V2)
correlation <- cor(quantile_normalized_df$logRe1, quantile_normalized_df$logRe2, method = "pearson")

# Create the plot with ggplot2
colors <- c(00,00, rev(viridis(200, option = "plasma"))[2:150])
values <- c(0, seq(0.01, 1, length.out = 200))


p <- ggplot(quantile_normalized_df, aes(x = logRe1, y = logRe2)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colors = colors, values = values) +  # Use the custom color scale
  theme_minimal() +  # Use a minimal theme
  labs(fill = "Density") +
  xlab("Replicate1") +
  ylab("Replicate2") +
  theme(text = element_text(size=40),
        axis.line = element_line(colour = "black")) + # Adding black axis lines
  #scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +  # Set x-axis breaks
  #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))+# Set y-axis breaks
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # Add diagonal dotted line
  # Label for the gradient scale

# Add correlation coefficient as text
p <- p + geom_text(aes(label = sprintf("rho = %.8f", correlation), x = Inf, y = Inf),
                   hjust = 1.1, vjust = 1.1, color = "black", size = 5)

ggsave(plot=p,paste0(OutFilePath,OutFileName,"_CorrelationPlot.pdf"), width=13, height=10)
