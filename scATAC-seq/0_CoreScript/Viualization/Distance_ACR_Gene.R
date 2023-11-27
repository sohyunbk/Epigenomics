library(GenomicRanges)
library(ggplot2)

GetDistance <- function(ACRfile,GeneAnn){
acr_data <- read.table(ACRfile, header = FALSE)
genes_data <- read.table(GeneAnn, header = FALSE)

# Convert to GRanges
acr_ranges <- GRanges(seqnames=acr_data$V1, 
                      ranges=IRanges(start=acr_data$V2, end=acr_data$V3))

genes_ranges <- GRanges(seqnames=genes_data$V1, 
                        ranges=IRanges(start=genes_data$V2, end=genes_data$V3),
                        gene_id=genes_data$V4)

# Find distance to nearest gene
nearest_genes <- distanceToNearest(acr_ranges, genes_ranges)
distances <- mcols(nearest_genes)$distance
nearest_genes_id <- genes_ranges[subjectHits(nearest_genes)]$gene_id

# Prepare the output
output <- data.frame(acr_data, nearest_genes_id, distances)
colnames(output) <- c("Chromosome", "Start", "End", "Closest_Gene_ID", "Distance")
return(output)
}

PieChart <- function(DistanceTable,OutputName){
  OverlapWithGene <- sum(DistanceTable$Distance == 0)
  NotOverlapWithGene <- sum(DistanceTable$Distance != 0)
  
  ## plot data
  df <- data.frame(
    category = c("Overlap with gene", "Intergenic"),
    value = c(OverlapWithGene, NotOverlapWithGene)
  )
  df$ratio <- df$value / sum(df$value)
  df$label_position <- cumsum(df$value) - df$value / 2
  df$label <- sprintf("%d (%.2f%%)",  df$value, df$ratio*100)
  
  ggplot(df, aes(x = "", y = value, fill = category)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y") +
    geom_text(aes(y = label_position, label = label), color = "white", size = 4) + 
    scale_fill_manual(values = c("orange", "purple")) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(OutputName, width=5, height=5)
}  
  
# Print the output
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn")
A619Higher <- "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC_FDR.0.01_A619Higher.Bed"
Bif3Higher <- "/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC_FDR.0.01_Bif3Higher.Bed"
AnnBed <- "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed"

A619Higher_GeneDistance <- GetDistance(A619Higher,AnnBed)
Bif3Higher_GeneDistance <- GetDistance(Bif3Higher,AnnBed)

PieChart(A619Higher_GeneDistance,"IM_OC_FDR.0.01_A619Higher.GeneDistancePieChart.pdf")
PieChart(Bif3Higher_GeneDistance,"IM_OC_FDR.0.01_Bif3Higher.GeneDistancePieChart.pdf")


NonZeroTable <- rbind(
  data.frame(Distance = A619Higher_GeneDistance$Distance[A619Higher_GeneDistance$Distance != 0], 
             dACR = "A619Higher"),
  data.frame(Distance = Bif3Higher_GeneDistance$Distance[Bif3Higher_GeneDistance$Distance != 0], 
             dACR = "Bif3Higher")
)
dim(NonZeroTable)
## Density plot

upper_limit <- 150000

ggplot(NonZeroTable, aes(x = Distance, fill = dACR)) +
  geom_density(alpha = 0.4) +
  coord_cartesian(xlim = c(0, upper_limit)) + # Set the x-axis limits
  theme_minimal() +
  scale_fill_manual(values = c("A619Higher" = "blue", "Bif3Higher" = "red")) + 
  labs(title = "Density Plot of Distance (Truncated)",
       x = "Distance",
       y = "Density")

ggsave("IM_OC_FDR.0.01_DensityPlot.pdf", width=7, height=5)

