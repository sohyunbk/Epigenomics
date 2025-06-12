##For quick I ran it through my mac

library(ggplot2)
Bif3Higher <- data.frame(Motif=c("ATCAT/ATGAT",
                                 "CGACGGCGAC/GTCGCCGTCG",
                                 "GCTAGC/GCTAGC",
                                 "AATWATT/AATWATT",
                                 "CATG/CATG"),
                         Ratio=c(37.7,
                                 14.9,
                                 21.9,
                                 49.4,
                                 54.7))
setwd("/Users/sohyun/Documents/2.SingleCellATAC/13.Data")
Bif3Higher$Motif <- reorder(Bif3Higher$Motif, Bif3Higher$Ratio)

ggplot(Bif3Higher, aes(y = Motif, x = Ratio)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Ratio", y = "Motif") +
  ggtitle("") +
  theme_minimal()+
  xlab("% of ACRs with the motif")+
  coord_cartesian(xlim = c(0, 60))

ggsave("Bif3HigherdACR_Motif_bar_plot.pdf", width = 4, height = 3)

######

A619Higher <- data.frame(Motif=c("ATCAT/ATGAT",
                                 "CGACGGCGAC/GTCGCCGTCG",
                                 "GCTAGC/GCTAGC",
                                 "AATWATT/AATWATT",
                                 "CATG/CATG"),
                         Ratio=c(37.7,
                                 14.9,
                                 21.9,
                                 49.4,
                                 54.7))
A619Higher$Motif <- reorder(A619Higher$Motif, A619Higher$Ratio)

ggplot(A619Higher, aes(y = Motif, x = Ratio)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Ratio", y = "Motif") +
  ggtitle("") +
  theme_minimal()+
  xlab("% of ACRs with the motif")+
  coord_cartesian(xlim = c(0, 60))

ggsave("A619HigherdACR_Motif_bar_plot.pdf", width = 4, height = 4)

