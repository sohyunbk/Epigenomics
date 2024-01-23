WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"
A619_meta <- paste0(WD,"/A619/Ref_AnnV4_metadata.txt")
loaded_A619_meta <- read.table(A619_meta)
bif3_meta <- paste0(WD,"/Bif3/Bif3_AnnV4_metadata.txt")
loaded_bif3_meta <- read.table(bif3_meta)
Slot <- "Ann_v4"
CluterTable <- rbind(loaded_A619_meta,loaded_bif3_meta)
setwd(WD)
head(CluterTable)
###########################
### 1) Barplot by celltype 
###########################
library(ggplot2)
source("/home/sb14489/Epigenomics/scATAC-seq/0_Function/DrawFigures_QC_Annotation_forUMAP.R")

Plotlist <- list()
ggplotData <- list()
k = 1
CellOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt")
colorr

##### 1) Separate bar plot

for (i in CellOrder){
  
  SubCluster <- CluterTable[which(CluterTable[[Slot]]==i),]
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
#str(ggplotData[[1]])
pdf(paste0(Slot,"_Ratio_by_Cluster_Barplot.pdf"),
     width = 10, height = 7, units = 'in', res = 300)
ggplot2.multiplot(ggplotData[[1]], ggplotData[[2]], ggplotData[[3]], ggplotData[[4]],
                  ggplotData[[5]], ggplotData[[6]], ggplotData[[7]], ggplotData[[8]],
                  ggplotData[[9]], ggplotData[[10]], ggplotData[[11]], ggplotData[[12]],
                  ggplotData[[13]],ggplotData[[14]],
                  cols=4)
dev.off()


###########################
### 2) StackedBarplot by Replicates 
###########################
Plotdata<- data.frame()
for (i in CellOrder){
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
Plotdata$Celltype  <- factor(Plotdata$Celltype,levels=CellOrder)
levels(factor(Plotdata$Celltype))
named_colors <- setNames(colorr, levels(Plotdata$Celltype))
ggplot(Plotdata, aes(fill=Celltype, y=Ratio, x=library)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = named_colors) +
  geom_text(aes(label = paste(round(Ratio,1),"%")), position = position_stack(vjust =  0.5))+
  theme_minimal()

ggsave("Ann_v3_Ratio_StackedBarplot_A619Re1and2_Bif3Re1and2.pdf", width=12, height=5)

###########################
### 3) StackedBarplot by Sample 
###########################
###########################
Plotlist <- list()
ggplotData <- list()


for (i in CellOrder){
  
  SubCluster <- CluterTable[which(CluterTable[[Slot]]==i),]
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
for (i in CellOrder){
  temp <- data.frame(Plotlist[[i]],Celltype = i)
  Plotdata <- rbind(Plotdata,temp)
}

library(plyr)
Plotdata <- ddply(Plotdata, .(SampleName),
                  transform, pos = cumsum(Ratio) - (0.5 * Ratio))

Plotdata$Ratio_Round <- round(Plotdata$Ratio, digits = 1)
Plotdata$Ratio_Round <- paste0(Plotdata$Ratio_Round," (",Plotdata$Fre,")")
tail(Plotdata)
Plotdata$Celltype  <- factor(Plotdata$Celltype,levels=CellOrder)

ggplot(Plotdata, aes(fill=Celltype, y=Ratio, x=SampleName)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colorr) +
  geom_text(aes(label = Ratio_Round), position = position_stack(vjust =  0.5))+
  theme_minimal()

ggsave(paste0(Slot,"_Ratio_StackedBarplot_BySample.pdf"), width=8, height=5)
###########################
### 4) Table save
###########################
SavedTable <- data.frame()

for (i in CellOrder){
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


write.table(SavedTable, row.names = TRUE,col.names=T,file=paste0(Slot,"CellNumber_A619_Bif3_Re12.txt"), sep="\t")
