setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn")
FileEnd <- ".EdgeRResult_PseudoReplicate_withPromoterRegion.txt"
Celltype <- 
  c("L1","L1atFloralMeristem",
    "FloralMeristem_SuppressedBract",
    "IM-OC","SPM-base_SM-base","IM_SPM_SM",
     "ProcambialMeristem_ProtoXylem_MetaXylem",
    "PhloemPrecursor", 
    "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma",
     "XylemParenchyma_PithParenchyma",
     "BundleSheath_VascularSchrenchyma",
    "CalloseRelated","G2_M")

Sig <- c()
Intergenic_ACR <- c()
Higher_Bif3 <- c()
Higher_WT <- c()
CutOffFDR <- 0.01
for (ct in CellTypeOrder){
  Data <- read.table(paste0(ct,FileEnd),header=TRUE)
  Sig <- c(Sig,sum(Data$FDR<CutOffFDR))
  Intergenic_ACR <- c(Intergenic_ACR,nrow(Data))
  Higher_Bif3 <- c(Higher_Bif3,sum(Data$FDR<CutOffFDR &Data$logFC >0))
  Higher_WT <- c(Higher_WT,sum(Data$FDR<CutOffFDR & Data$logFC <0))
}

dACRInfo <- data.frame(Celltype,Sig,Intergenic_ACR)
Bif3Higher  <- data.frame(Celltype, Sig=Higher_Bif3,Intergenic_ACR)
WTHigher  <- data.frame(Celltype,Sig=Higher_WT,Intergenic_ACR)
dACRInfo$dACRRatio <- (dACRInfo$Sig / dACRInfo$Intergenic_ACR)*100
Bif3Higher$dACRRatio <- (Bif3Higher$Sig / Bif3Higher$Intergenic_ACR)*100
WTHigher$dACRRatio <- (WTHigher$Sig / WTHigher$Intergenic_ACR)*100

dACRInfo$Color <- "dACRTotal"
Bif3Higher$Color <- "Bif3Higher"
WTHigher$Color <- "WTHigher"
head(dACRInfo)
head(Bif3Higher)

FigureTable <- rbind(dACRInfo,Bif3Higher,WTHigher)

library(ggplot2)
custom_colors <- c("#4d0505","#092691", "#ad8c09") 
FigureTable$Color <- factor(FigureTable$Color,levels=c("dACRTotal","Bif3Higher","WTHigher"))
ggplot(FigureTable, aes(x = Sig, y = Celltype, size = dACRRatio)) +
  geom_point(aes(colour = factor(Color)),alpha=0.4)+  # New geom_point for Higher_Bif3Ratio
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = custom_colors) + 
  theme_minimal() +
  labs(x = "The number of dACR", y = "Cell type", title = " ", size = "dACR Ratio") +
  theme(axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),  # Adjust legend text size if needed
        legend.title = element_text(size = 8), # Adjust legend title size if needed
        axis.title.y = element_text(size = 8), # Adjust y-axis title size if needed
        axis.title.x = element_text(size = 8), # Adjust x-axis title size if needed
        title = element_text(size = 12), # Adjust title size if needed
        plot.title = element_text(size = 14), # Adjust plot title size if needed
        strip.text = element_text(size = 8), # Adjust facet strip text size if needed
        axis.ticks = element_blank(), # Remove axis tick marks
        axis.line = element_line(colour = "black"), # Add axis lines
        axis.line.x = element_line(), # Customize x-axis line separately
        axis.line.y = element_line() # Customize y-axis line separately
  )+labs(color = " ")

ggsave("dACRNumber_FDR0.01.pdf"
       , width=7, height=3)


