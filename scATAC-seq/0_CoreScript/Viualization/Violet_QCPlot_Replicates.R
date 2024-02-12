#### Ann V4 version!
library(ggplot2)
library(tidyr)

read_and_prepare_data <- function(file_path) {
  # Read the data from the file
  meta_data <- read.table(file_path)
  
  # Select specific columns and create a new data frame
  data_frame <- data.frame(
    library = meta_data$library, 
    FRiP = meta_data$FRiP,
    pPtMt = meta_data$pPtMt, 
    pTSS = meta_data$pTSS,
    log10nSites = meta_data$log10nSites,
    Ann = meta_data$Ann_v4
  )
  
  return(data_frame)
}

# Using the function to read and prepare the data
A619Meta_path <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt"
Bif3Meta_path <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV4_metadata.txt"

A619Meta <- read_and_prepare_data(A619Meta_path)
Bif3Meta <- read_and_prepare_data(Bif3Meta_path)

# Combining the data frames
PlotData <- rbind(A619Meta, Bif3Meta)

head(PlotData)
tail(PlotData)

CellOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt")
# Reshape the data to a long format
PlotDataLong <- pivot_longer(PlotData, cols = c(FRiP, pPtMt, pTSS, log10nSites), 
                             names_to = "Metric", values_to = "Value")
PlotDataLong<- as.data.frame(PlotDataLong)
PlotDataLong$Ann <- factor(PlotDataLong$Ann,levels=c(CellOrder))
ggplot(PlotDataLong, aes(x = Ann, y = Value, fill = library)) +
  geom_violin(trim = FALSE, color = NA) + 
  facet_wrap(~Metric, scales = "free_y", nrow = 1) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("A619_Re1" = "#71f562", "A619_Re2" = "#d5eb6c", "bif3_Re1" ="#6cebeb" , "bif3_Re2" = "#7392fa")) +
  labs(fill = "Library")

ggsave("/scratch/sb14489/3.scATAC/2.Maize_ear/9.CheckQC/1.Tn5Numb_TSS_FRiP_OrgRatio_byCellTypes/A619_Bif3_Replicates_Tn5_FRiP_TSS_Org.pdf"
  , width=23, height=6)	


