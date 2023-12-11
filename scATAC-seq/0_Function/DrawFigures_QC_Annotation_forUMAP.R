library(ggplot2)
library(stringr)
library("RColorBrewer")
library(gridExtra)

### Getting metafile as input.
DrawUMAP_Ann_QC <- function(PreAnnMeta,Meta, Slot, CellOrder, Re1, Re2,OutfilePathName){
  print("Replicates name should be sampleID slot")
ColorForPreAnn <- c( "#f58c8c", "#FFB380", "#FFECB3", "#a2d179", "#B3FFB3","#44c78e", "#B3FFE5", 
                     "#8dccca", "#48a2f0", "#69a0cf", "#b2bcf7", "#D1B3FF", "#ECB3FF",
                     "#a088f7")
PreAnn <- ggplot(PreAnnMeta, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
    geom_point(size=0.02) +
    scale_color_manual(values=ColorForPreAnn)+theme_minimal()+
    guides(colour = guide_legend(override.aes = list(size=7)))+
    labs(title = paste0("Pre-Annotation\n CellNumber: ",nrow(PreClusterMeta)),
         x = "UMAP1",
         y = "UMAP2")
    
Meta$Ann <- Meta[[Slot]]
Meta$Ann <- factor(Meta$Ann,levels=(CellOrder))
colorr <- c("#4F96C4","#84f5d9","#0bd43d","#d62744","#FDA33F","#060878","#62a888",
            "#876b58","#800000", "#800075","#e8cf4f","#adafde","#DE9A89","#5703ff",
            "#deadce","#fc53b6")

All <- ggplot(Meta, aes(x=umap1, y=umap2, color=factor(Ann))) +
  geom_point(size=0.02) +
  scale_color_manual(values=colorr)+theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  labs(title = paste0("Re1+R2 \n CellNumber: ",nrow(Meta)),
       x = "UMAP1",
       y = "UMAP2")

ClustersTable_Re1 <- subset(Meta, sampleID == Re1)
Re1_plot <- ggplot(ClustersTable_Re1, aes(x=umap1, y=umap2, color=factor(Ann))) +
  geom_point(size=0.02, color="blue") +
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(title = paste0("Re1 : ",nrow(ClustersTable_Re1)))
ClustersTable_Re2 <- subset(Meta, sampleID == Re2)
Re2_plot <- ggplot(ClustersTable_Re2, aes(x=umap1, y=umap2, color=factor(Ann))) +
  geom_point(size=0.02, color="red") +
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(title = paste0("Re2 : ",nrow(ClustersTable_Re2)))

## Add some plots to see the cell quality
### * Tn5 log
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(Meta$log10nSites),
                                      max(Meta$log10nSites)))
Q_Tn5 <- ggplot(Meta, aes(x=umap1, y=umap2,
                                color=log10nSites)) +
  geom_point(size=0.02) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "LogTn5")+sc
## * Doublets
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(Meta$doubletscore),
                                      max(Meta$doubletscore)))
Q_doubletscore <-ggplot(Meta, aes(x=umap1, y=umap2,
                                        color=doubletscore)) +
  geom_point(size=0.02) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "Doublet score")+sc
## Tss ratio
head(Meta)
Meta$rTSS <- Meta$tss/Meta$total
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(Meta$rTSS),
                                      max(Meta$rTSS)))
Q_rTSS <-ggplot(Meta, aes(x=umap1, y=umap2,
                                color=rTSS)) +
  geom_point(size=0.02) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "TSS ratio")+sc

## Tss ratio
head(Meta)
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(Meta$FRiP),
                                      max(Meta$FRiP)))
Q_FRiP <-ggplot(Meta, aes(x=umap1, y=umap2,
                          color=FRiP)) +
  geom_point(size=0.02) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "FRiP ratio")+sc

pdf(paste0(OutfilePathName,"_AnnQCPlots.pdf"), width=44, height=5)
grid.arrange(PreAnn,All, Re1_plot, Re2_plot,
             Q_Tn5, Q_doubletscore,Q_rTSS,Q_FRiP,
             ncol=8, widths=c(1.5,2.1,1,1,1.2,1.2,1.2,1.2))
dev.off()
}