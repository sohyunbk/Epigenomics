library(dplyr)
library(GenomicRanges)

#### ******* Test dataset
SNPFile_test <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs.v2.txt", header=TRUE)
unique_SNPFile_test <- unique(SNPFile_test)
print(nrow(SNPFile_test))
print(nrow(unique_SNPFile_test))
head(SNPFile_test)
element_counts <- table(unique_SNPFile_test$acrID)
head(element_counts)
SNPFile_test[SNPFile_test$snpID == "chr9_99008396",]
head(unique(SNPFile_test$acrID))
length(unique(SNPFile_test$acrID))
## 1) Check the overlap as it's extended with 1000bp
split_pos <- strsplit(unique(SNPFile_test$acrID), "_")
chr <- sapply(split_pos, function(x) x[1])
start <- as.integer(sapply(split_pos, function(x) x[2])) -500
end <- as.integer(sapply(split_pos, function(x) x[3])) +500
gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

overlaps <- findOverlaps(gr)
summary(overlaps)
overlap_pairs <- as.data.frame(overlaps)
overlap_counts <- countOverlaps(gr)
NonOverlappedACRs <- unique(SNPFile_test$acrID)[overlap_counts==1]
head(NonOverlappedACRs)

NonOverlappedACRs_split <- strsplit(NonOverlappedACRs, "_")
chr <- sapply(NonOverlappedACRs_split, function(x) x[1])
start <- as.integer(sapply(NonOverlappedACRs_split, function(x) x[2])) -250
end <- as.integer(sapply(NonOverlappedACRs_split, function(x) x[3])) +250
BedFile <- data.frame(chr,start,end)
head(BedFile)
BedFile$Strand <- "+"
BedFile$CellTypeTemp <- "1"
dim(BedFile)
write.table(BedFile, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs.v2.curated1000bp.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=F)

dim(SNPFile_test)
SNPFile_test_Filtered <- SNPFile_test[SNPFile_test$acrID %in% NonOverlappedACRs,]
dim(SNPFile_test_Filtered)
UniqueACR_Test <- SNPFile_test_Filtered %>%
  group_by(acrID) %>%
  sample_n(size = 1) %>%
  ungroup()
dim(UniqueACR_Test)

write.table(UniqueACR_Test, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs.v2.curated.txt", sep = "\t", row.names = FALSE, quote=F)

#### ****** Control File
### Check common ACRs
SNPFile <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs.v2.txt", header=TRUE)
unique_SNPFile <- unique(SNPFile)
print(nrow(SNPFile))
print(nrow(unique_SNPFile))
head(unique_SNPFile)
## 1) Check the overlap as it's extended with 1000bp
split_pos <- strsplit(unique(SNPFile$acrID), "_")
chr <- sapply(split_pos, function(x) x[1])
start <- as.integer(sapply(split_pos, function(x) x[2])) -500
end <- as.integer(sapply(split_pos, function(x) x[3])) +500
gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
overlaps <- findOverlaps(gr)
summary(overlaps)
overlap_pairs <- as.data.frame(overlaps)
overlap_counts <- countOverlaps(gr)
NonOverlappedACRs <- unique(SNPFile_test$acrID)[overlap_counts==1]
head(NonOverlappedACRs)
dim(SNPFile_test)
SNPFile_Filtered <- SNPFile[SNPFile$acrID %in% NonOverlappedACRs,]
dim(SNPFile_Filtered)

UniqueACR <- SNPFile_Filtered %>%
  group_by(acrID) %>%
  sample_n(size = 1) %>%
  ungroup()

write.table(UniqueACR, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs.v2.curated.txt", sep = "\t", row.names = FALSE, quote=F)
