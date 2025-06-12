dACRInfo <- data.frame(
  "Celltype" = c(
    "BundleSheath or VascularSchrenchyma",
    "CalloseRelated",
    "FloralMeristem or SuppressedBract",
    "G2&M",
    "IM-OC",
    "IM or SPM or SM",
    "L1",
    "L1atFloralMeristem",
    "PhloemPrecursor",
    "Proto/MetaXylem or ProcambialMeristem",
    "Proto/MetaPhloem or CC or PP",
    "SPM-base or SM-base",
    "XylemParenchyma or PithParenchyma"
  ),
  "Sig" = c(380, 29, 58, 80, 3632, 296, 46, 187, 113, 222, 10, 446, 433),
  "Intergenic_ACR" = c(48136, 39022, 36261, 35570, 37306,
                       51758, 38363, 42950, 48565, 39351, 26873, 38290, 45334)
)

dACRInfo
dACRInfo$dACRRatio <- (dACRInfo$Sig / dACRInfo$Intergenic_ACR)*100


library(ggplot2)

ggplot(dACRInfo, aes(x = Sig, y = Celltype, size = dACRRatio)) +
  geom_point(color = "brown") +
  scale_size_continuous(range = c(1, 5)) +
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
  )

ggsave("/Users/sohyun/Documents/2.SingleCellATAC/dACRNumber.pdf"
       , width=5.3, height=3)


ggplot(dACRInfo, aes(x = reorder(`Celltype`, Sig), y = Sig)) +
  geom_point(size = 3, color = "blue") +
  geom_segment(aes(xend = reorder(`Celltype`, Sig), yend = 0), color = "gray50") +
  labs(x = "Cell type", y = "Number of Sig", title = "Number of Significant Values by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


