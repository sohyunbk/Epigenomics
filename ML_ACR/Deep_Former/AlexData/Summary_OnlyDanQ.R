
DanQ_MappableRegion <- read.table("/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative/test_performance.txt",
                            header=TRUE)
library(ggplot2)
setwd("/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative")
head(DanQ_MappableRegion)

ggplot(data=DanQ_MappableRegion, aes(x=class, y=roc_auc)) +
  geom_bar(stat="identity")+
  theme_minimal()+ 
  #scale_fill_manual(values=c("#E69F00"))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, 
                                   size = 10, hjust = 1))+
  xlab(" ") 

ggsave("MappableGenome.pdf", width=15, height=5)
