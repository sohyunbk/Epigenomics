setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_SplitPeaks")
Distance <- read.table("Distance.txt")
library(ggplot2)
ggplot(Distance, aes(x=V2)) + 
  geom_density()+
  scale_x_continuous(n.breaks=20)+
  xlab("Distance")

ggsave("DensityPlot.pdf", width=12, height=10)


setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_SplitPeaks")
Distance <- read.table("PeakLength")
head(Distance)
library(ggplot2)
ggplot(Distance, aes(x=V1)) + 
  geom_density()+
  scale_x_continuous(n.breaks=20)+
  xlab("Length")

ggsave("DensityPlot_Length.pdf", width=12, height=10)


mean(Distance$V1)
median(Distance$V1)
