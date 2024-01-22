A619_meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
loaded_A619_meta <- read.table(A619_meta)
bif3_meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt"
loaded_bif3_meta <- read.table(bif3_meta)

CluterTable <- rbind(loaded_A619_meta,loaded_bif3_meta)
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3")
head(CluterTable)
###########################
### 1) Barplot by celltype 
###########################
library(ggplot2)
Plotlist <- list()
ggplotData <- list()
i <- levels(factor(CluterTable$Ann_v3))[1]
k = 1

colorr <-c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
           "#876b58","#800000", "#800075","#e8cf4f","#0bd43d","#fc53b6",
           "#deadce","#adafde","#5703ff")

for (i in levels(factor(CluterTable$Ann_v3))){
  
  SubCluster <- CluterTable[which(CluterTable$Ann_v3==i),]
  #head(SubCluster)
  Plotlist[[i]] <- data.frame()
  for (j in levels(factor(SubCluster$library))){
    Count <- nrow(SubCluster[which(SubCluster$library==j),])
    Temp <- data.frame(library=j,Fre=Count,
                       Ratio=(Count/table(CluterTable$library)[[j]])*100)
    #rownames(Temp) <- NULL
    Plotlist[[i]] <- rbind(Plotlist[[i]],Temp)
  }
  Plotlist[[i]]$library <- factor(Plotlist[[i]]$library ,levels=c("A619_Re1", "A619_Re2", "bif3_Re1","bif3_Re2"
  ))
  ggplotData[[i]] <- ggplot(data=Plotlist[[i]], aes(x=library, y=Ratio)) +theme_minimal()+
    geom_bar(stat="identity",fill=colorr[k])+coord_flip()+
    xlab("") +
    ylab("Ratio(%)")+ggtitle(i)+ theme(plot.title = element_text(size = 10))
  k=k+1
}

library(easyGgplot2)
getwd()
head(SubCluster)
str(ggplotData[[1]])
tiff("Ann_v3_Ratio_by_Cluster_Barplot.tiff", width = 10, height = 7, units = 'in', res = 300)
ggplot2.multiplot(ggplotData[[1]], ggplotData[[2]], ggplotData[[3]], ggplotData[[4]],
                  ggplotData[[5]], ggplotData[[6]], ggplotData[[7]], ggplotData[[8]],
                  ggplotData[[9]], ggplotData[[10]], ggplotData[[11]], ggplotData[[12]],
                  ggplotData[[13]],
                  cols=4)
dev.off()


###########################
### 2) StackedBarplot by Replicates 
###########################
Plotdata<- data.frame()
for (i in levels(factor(CluterTable$Ann_v3))){
  temp <- data.frame(Plotlist[[i]],Celltype = i)
  Plotdata <- rbind(Plotdata,temp)
}

head(Plotdata)
Plotdata$Celltype <- factor(Plotdata$Celltype)
Plotdata$library <- factor(Plotdata$library)
library(plyr)
Plotdata <- ddply(Plotdata, .(library),
                  transform, pos = cumsum(Ratio) - (0.5 * Ratio))

Plotdata$Ratio_Round <- round(Plotdata$Ratio, digits = 1)
tail(Plotdata)

ggplot(Plotdata, aes(fill=Celltype, y=Ratio, x=library)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colorr[1:13]) +
  geom_text(aes(label = paste(round(Ratio,1),"%")), position = position_stack(vjust =  0.5))+
  theme_minimal()

ggsave("Ann_v3_Ratio_StackedBarplot_A619Re1and2_Bif3Re1and2.pdf", width=12, height=5)

###########################
### 3) StackedBarplot by Sample 
###########################
###########################
Plotlist <- list()
ggplotData <- list()
i <- levels(factor(CluterTable$Ann_v3))[1]
k = 1

colorr <-c("#4F96C4","#84f5d9","#DE9A89","#FDA33F","#060878","#d62744","#62a888",
           "#876b58","#800000", "#800075","#e8cf4f","#0bd43d","#fc53b6",
           "#deadce","#adafde","#5703ff")
head(CluterTable)
for (i in levels(factor(CluterTable$Ann_v3))){
  
  SubCluster <- CluterTable[which(CluterTable$Ann_v3==i),]
  #head(SubCluster)
  Plotlist[[i]] <- data.frame()
  for (j in levels(factor(SubCluster$SampleName))){
    Count <- nrow(SubCluster[which(SubCluster$SampleName==j),])
    Temp <- data.frame(SampleName=j,Fre=Count,
                       Ratio=(Count/table(CluterTable$SampleName)[[j]])*100)
    #rownames(Temp) <- NULL
    Plotlist[[i]] <- rbind(Plotlist[[i]],Temp)
  }
  Plotlist[[i]]$SampleName <- factor(Plotlist[[i]]$SampleName ,levels=c("A619", "bif3"
  ))
  ggplotData[[i]] <- ggplot(data=Plotlist[[i]], aes(x=SampleName, y=Ratio)) +theme_minimal()+
    geom_bar(stat="identity",fill=colorr[k])+coord_flip()+
    xlab("") +
    ylab("Ratio(%)")+ggtitle(i)+ theme(plot.title = element_text(size = 10))
  k=k+1
}

Plotdata<- data.frame()
for (i in levels(factor(CluterTable$Ann_v3))){
  temp <- data.frame(Plotlist[[i]],Celltype = i)
  Plotdata <- rbind(Plotdata,temp)
}

library(plyr)
Plotdata <- ddply(Plotdata, .(SampleName),
                  transform, pos = cumsum(Ratio) - (0.5 * Ratio))

Plotdata$Ratio_Round <- round(Plotdata$Ratio, digits = 1)
Plotdata$Ratio_Round <- paste0(Plotdata$Ratio_Round," (",Plotdata$Fre,")")
tail(Plotdata)

ggplot(Plotdata, aes(fill=Celltype, y=Ratio, x=SampleName)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colorr) +
  geom_text(aes(label = Ratio_Round), position = position_stack(vjust =  0.5))+
  theme_minimal()

ggsave("Ann_v3_Ratio_StackedBarplot_BySample.pdf", width=8, height=5)
###########################
### 4) Table save
###########################
SavedTable <- data.frame()

for (i in levels(factor(CluterTable$Ann_v3))){
  # Order = A619_Re1, A619_Re2, bif3_Re1, bif3_Re2
  temp <- data.frame(Plotlist[[i]],Celltype = i)
  tTable <- data.frame(t(c(temp$Fre)))
  colnames(tTable) <- c("A619_Re1","A619_Re2","Bif3_Re1","Bif3_Re2")
  rownames(tTable) <- i
  SavedTable <- rbind(SavedTable,tTable)
}

ColSum <- t(data.frame(colSums(SavedTable)))
rownames(ColSum) <- c("Sum")
SavedTable <- rbind(SavedTable, ColSum)

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/13.CheckQC")

write.table(SavedTable, row.names = TRUE,col.names=T,file="CellNumber_A619_Bif3_Re12.txt", sep="\t")
