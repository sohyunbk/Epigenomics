##For quick I ran it through my mac

######## This is the previous version, and it's from STEME #####
#### However,  the number of motif occurrence from STREME and fimo are different... I am following fimo

setwd("/Users/sohyun/Documents/2.SingleCellATAC/")

library(ggplot2)
Bif3Higher <- data.frame(Motif=c("CAATAATGC"),
                         Ratio=c(23.4))

Bif3Higher$Fill <- "Bif3Higher"

######


A619Higher <- data.frame(Motif=c("GCACAGCAGCR",
                                 "GCAGCATGC",
                                 "YAGAGAGAGA",
                                 "CGCGCCGCGCC",
                                 "GCTAGCTAGC"
                                 ),
                         Ratio=c(29.4,
                                 24.8,
                                 17.9,
                                 5.7,
                                 11.7 ))



A619Higher <- data.frame(Motif=c("GCACAGCAGCR",
                                 "GCAGCATGC",
                                 "YAGAGAGAGA",
                                 "CGCGCCGCGCC",
                                 "GCTAGCTAGC",
                                 "CAGTGGCA",
                                 "ATGCATGCAAT",
                                 "GCAGGCAGGC"),
                         Ratio=c(29.4,
                                 24.8,
                                 17.9,
                                 5.7,
                                 11.7,
                                 26.5,
                                 18.7,
                                 19.3 ))

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

#################################################################
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4")
Fimo_TAAT <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/fimo_out_2/fimo.tsv",header=TRUE)


