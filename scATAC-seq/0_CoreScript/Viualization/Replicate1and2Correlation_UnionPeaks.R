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
  make_option(c("--Re1_BulkPeak"), type="character",
              help="Re1_BulkPeak", metavar="character"),
  make_option(c("--Re2_BulkPeak"), type="character",
              help="Re2_BulkPeak"),
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


#Replicate1_BulkPeak <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619/macs2_temp/bulk_peaks_peaks.narrowPeak"
#Replicate2_BulkPeak <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619_2/macs2_temp/bulk_peaks_peaks.narrowPeak"
#Replicate1_AllReads <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Unique.bed"
#Replicate2_AllReads <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_2_Unique.bed"
Replicate1_BulkPeak <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/3_bif3/macs2_temp/bulk_peaks_peaks.narrowPeak"
Replicate2_BulkPeak <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/3_bif3_2/macs2_temp/bulk_peaks_peaks.narrowPeak"
Replicate1_AllReads <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Unique.bed"
Replicate2_AllReads <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_2_Unique.bed"

Replicate1_BulkPeak <- opt$Re1_pBulkPeak
Replicate2_BulkPeak <- opt$Re2_BulkPeak
Replicate1_AllReads  <- opt$Re1_AllReads
Replicate2_AllReads <- opt$Re2_AllReads
OutFilePath <- opt$OutPath
OutFileName <- opt$OutFileName


## Function1: Read Summit
ReadPeak <- function(PeakInput){

PeakFile <- read.table(PeakInput)
#head(PeakFile)
PeakFile_chr <- PeakFile[grepl("^chr", PeakFile$V1), ]
PeakFile_Grange <-  GRanges(seqnames = PeakFile_chr$V1,
                       ranges = IRanges(start = PeakFile_chr$V2,
                                        end = PeakFile_chr$V3
                                    ))

return(PeakFile_Grange)
}

Peak1 <- ReadPeak(Replicate1_BulkPeak)
Peak2 <- ReadPeak(Replicate2_BulkPeak)
CommonPeak_GRange_unique <- reduce(c(Peak1, Peak2)) ## Union

Overlap_Tn5_CommonPeak <- function(Tn5File, CommonPeak_GRange_unique) {
  print("Read Tn5 file")
  AllTn5 <- read.table(Tn5File, stringsAsFactors = FALSE)
  print("Done: Read Tn5 file")
  
  Tn5_GRange <- GRanges(seqnames = AllTn5$V1,
                        ranges = IRanges(start = AllTn5$V2, end = AllTn5$V3),
                        names = AllTn5$V4)
  
  print('Getting overlap between common peak and Tn5')
  overlaps_final <- findOverlaps(Tn5_GRange, CommonPeak_GRange_unique)
  print('Done: getting overlap between common peak and Tn5')
  
  # Efficiently count overlaps
  overlap_counts <- table(subjectHits(overlaps_final))
  
  # Create a data frame for the results
  result <- data.frame(seqnames = seqnames(CommonPeak_GRange_unique)[as.integer(names(overlap_counts))],
                       start = start(CommonPeak_GRange_unique)[as.integer(names(overlap_counts))],
                       end = end(CommonPeak_GRange_unique)[as.integer(names(overlap_counts))],
                       count = overlap_counts)
  
  return(result)
}



overlap_counts_Re1 <- Overlap_Tn5_CommonPeak(Replicate1_AllReads,CommonPeak_GRange_unique)
overlap_counts_Re2 <- Overlap_Tn5_CommonPeak(Replicate2_AllReads,CommonPeak_GRange_unique)
overlap_counts_Re2$`count.Freq_Re2` <- overlap_counts_Re2$count.Freq
overlap_counts_Re1$`count.Freq_Re1` <- overlap_counts_Re1$count.Freq
overlap_counts_Re2 <- subset(overlap_counts_Re2, select = -c(count.Freq, count.Var1))
overlap_counts_Re1 <- subset(overlap_counts_Re1, select = -c(count.Freq, count.Var1))
combined_table <- merge(overlap_counts_Re1, overlap_counts_Re2, by = c("seqnames", "start", "end"), all = TRUE)
combined_table$count.Freq_Re2[is.na(combined_table$count.Freq_Re2)] <- 0
combined_table$count.Freq_Re1[is.na(combined_table$count.Freq_Re1)] <- 0


CorrTable <- data.frame(peak=paste(combined_table$seqnames,
                                   combined_table$start,
                                   combined_table$end,
                                   sep="_"),
                        Re1_count=combined_table$count.Freq_Re1,
                        Re2_count=combined_table$count.Freq_Re2)

rownames(CorrTable) <- CorrTable$peak
CorrTable <- CorrTable[,-1]
# CPM Normalization
#quantile_normalized <- normalize.quantiles(as.matrix(CorrTable))
#quantile_normalized_df <- as.data.frame(quantile_normalized)
#quantile_normalized_df$logRe1 <- log10(quantile_normalized_df$V1)
#quantile_normalized_df$logRe2 <- log10(quantile_normalized_df$V2)
calculateCPM <- function(countData) {
  totalCounts <- colSums(countData)
  cpm <- sweep(countData, 2, totalCounts, "/") * 1e6
  return(cpm)
}
CPMTable <- calculateCPM(CorrTable)
correlation <- cor(CPMTable$Re1_count, CPMTable$Re2_count, method = "pearson")
correlation <- cor(CPMTable$Re1_count, CPMTable$Re2_count, method = "spearman")


# Create the plot with ggplot2
colors <- c(00,00, rev(viridis(200, option = "viridis"))[2:150])
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
