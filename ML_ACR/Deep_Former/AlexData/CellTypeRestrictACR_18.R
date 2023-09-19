NewACR_byCT <- read.table("/scratch/sb14489/8.ML_ACR/1.InputBed/celltype_ACRs.500bp.txt")
head(NewACR_byCT)
rowSums(NewACR_byCT)


## Check Histogram graph
library(stringr)
library(ggplot2)
library(tidyverse) 

head(as.numeric(rowSums(NewACR_byCT)))
ggplot() + aes(as.numeric(rowSums(NewACR_byCT)))+ geom_histogram(bins=31)
ggsave("/scratch/sb14489/8.ML_ACR/1.InputBed/FrequencyPlot_EditedFile.pdf", width=13, height=10)	

length(rownames(NewACR_byCT))
CellTypeRestrict <- rownames(NewACR_byCT)[rowSums(NewACR_byCT) < 33]
CellTypeRestrict[1]

## Write Bed File
NewACR_byCT$Pos <- rownames(NewACR_byCT)
head(NewACR_byCT)
NewACR_byCT_Indexed <- pivot_longer(NewACR_byCT, - Pos)
head(NewACR_byCT_Indexed)

NewACR_byCT_Indexed_OnlyOne <- NewACR_byCT_Indexed[NewACR_byCT_Indexed$value>0,]
head(NewACR_byCT_Indexed_OnlyOne)
dim(NewACR_byCT_Indexed_OnlyOne)


### 1) All Cell type
NewACR_byCT_Indexed_OnlyOne
NewACR_byCT_Indexed_OnlyOne$chr <- as.character(lapply(strsplit(as.character(NewACR_byCT_Indexed_OnlyOne$Pos),
                                                                split="_"), "[", 1))
head(NewACR_byCT_Indexed_OnlyOne)
NewACR_byCT_Indexed_OnlyOne$start <- as.character(lapply(strsplit(as.character(NewACR_byCT_Indexed_OnlyOne$Pos),
                                                                split="_"), "[", 2))

NewACR_byCT_Indexed_OnlyOne$end <- as.character(lapply(strsplit(as.character(NewACR_byCT_Indexed_OnlyOne$Pos),
                                                                  split="_"), "[", 3))

head(NewACR_byCT_Indexed_OnlyOne)
BedFormat <- subset(NewACR_byCT_Indexed_OnlyOne,select=c("chr","start","end","name"))
BedFormat
dim(BedFormat)

write.table(BedFormat,"/scratch/sb14489/8.ML_ACR/1.InputBed/AllACRs_130539.500bp.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

### 
head(BedFormat)
head(NewACR_byCT_Indexed_OnlyOne)
NewACR_byCT_Indexed_OnlyOne_CellTypeRestrict <- NewACR_byCT_Indexed_OnlyOne[NewACR_byCT_Indexed_OnlyOne$Pos%in%CellTypeRestrict,]
NewACR_byCT_Indexed_OnlyOne_CellTypeRestrict
BedFormat_CTRestrict <- subset(NewACR_byCT_Indexed_OnlyOne_CellTypeRestrict,select=c("chr","start","end","name"))

write.table(BedFormat_CTRestrict,"/scratch/sb14489/8.ML_ACR/1.InputBed/celltype_RestrictACRs_118160.500bp.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

