library(ggplot2)



DrawPieChart <-function(WD,OutfileName,IntergenicName,GenicName){
  
IntergenicFile <- paste0(WD,IntergenicName)
GenicFile <-  paste0(WD,GenicName)
nIntergenic <- length(readLines(IntergenicFile))
nGenic <- length(readLines(GenicFile))
total <- nGenic + nIntergenic


# Create a data frame
data <- data.frame(
  category = c("Genic", "Intergenic"),
  count = c(nGenic, nIntergenic)
)

# Add a column for ratio
data$ratio <- data$count / total

# Add a column for labels (ratio and count)
data$labels <- paste0(round(data$ratio * 100), "% (", data$count, ")")

ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  geom_text(aes(label = labels, y = count / 2), color = "white") +
  theme_void() +
  theme(legend.title = element_blank()) +
  labs(fill = "Category")+
  scale_fill_manual(values = c("#9269c7", "#4dd6af")) 
  
ggsave(paste0(WD,OutfileName), width=4, height=5)
}

WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619/"
IntergenicName <- "A619.500bp_peaks_Intergenic.bed"
GenicName <- "A619.500bp_peaks_Genic.bed"
OutfileName <- "/A619_genic_intergenic_piechart.pdf"
DrawPieChart(WD,OutfileName,IntergenicName,GenicName)

WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/Bif3/"
IntergenicName <- "Bif3.500bp_peaks_Intergenic.bed"
GenicName <- "Bif3.500bp_peaks_Genic.bed"
OutfileName <- "/Bif3_genic_intergenic_piechart.pdf"
DrawPieChart(WD,OutfileName,IntergenicName,GenicName)
