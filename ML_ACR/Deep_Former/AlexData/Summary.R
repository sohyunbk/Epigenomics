
DeeperDeepSea_WholeGenome <- read.table("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/204bp_AllGenome/test_performance.txt",
                            header=TRUE)
DeeperDeepSea_MappableRegion <- read.table("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/204bp_MappableRegions/test_performance.txt",
                            header=TRUE)
DanQ_WholeGenome <- read.table("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome_DanQ/test_performance.txt",
                                        header=TRUE)
DanQ_MappableRegion <- read.table("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_MappableRegions_DanQ/test_performance.txt",
                                           header=TRUE)
head(DanQ_WholeGenome)

DeeperDeepSea_WholeGenome$Method <- "DeeperDeepSea"
DanQ_WholeGenome$Method <- "DanQ"

WholeGenome <- rbind(DeeperDeepSea_WholeGenome,DanQ_WholeGenome)
head(WholeGenome)

library(ggplot2)
setwd("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test")

ggplot(data=WholeGenome, aes(x=class, y=roc_auc, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+ 
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  
  theme(axis.text.x = element_text(angle = 60, vjust = 1, 
                                   size = 10, hjust = 1))+
  xlab(" ") 

ggsave("wholeGenome.pdf", width=15, height=5)


DeeperDeepSea_MappableRegion$Method <- "DeeperDeepSea"
DanQ_MappableRegion$Method <- "DanQ"

Mappable <- rbind(DeeperDeepSea_MappableRegion,DanQ_MappableRegion)
head(Mappable)

library(ggplot2)
setwd("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test")

ggplot(data=Mappable, aes(x=class, y=roc_auc, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+ 
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  
  theme(axis.text.x = element_text(angle = 60, vjust = 1, 
                                   size = 10, hjust = 1))+
  xlab(" ") 

ggsave("Mappable.pdf", width=15, height=5)
