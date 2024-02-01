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
dim(SNPFile_test)
SNPFile_test_Filtered <- SNPFile_test[SNPFile_test$acrID %in% NonOverlappedACRs,]
dim(SNPFile_test_Filtered)
UniqueACR_Test <- SNPFile_test_Filtered %>%
  group_by(acrID) %>%
  sample_n(size = 1) %>%
  ungroup()

write.table(UniqueACR_Test, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs.v2.curated.txt", sep = "\t", row.names = FALSE, quote=F)

#### 2) Make Input bed file! :)
### Check common ACRs
head(randomly_selected_control$acrID)
head(ACRs_pos1NonRe) 
head(randomly_selected_control$acrID)
length(intersect(randomly_selected_control$acrID, ACRs_pos_NonRe))
length(intersect(merged_SNPFile_test$acrID, ACRs_pos_NonRe))
SelectedACR_Control <- intersect(randomly_selected_control$acrID, ACRs_pos_NonRe)
SelectedACR_Test <- intersect(merged_SNPFile_test$acrID, ACRs_pos_NonRe)
sum(!grepl("^chr10", SelectedACR_Test))
SelectedACR_Control_SameNumberwithTest <- sample(SelectedACR_Control, size = length(intersect(merged_SNPFile_test$acrID, ACRs_pos_NonRe)))
head(SelectedACR_Control_SameNumberwithTest)

head(ACRs_18cells)
ACRs_18cells$V5 <- paste(ACRs_18cells$V1, ACRs_18cells$V2, ACRs_18cells$V3, sep="_")
bedfile_SelectedACR_Test <- ACRs_18cells[ACRs_18cells$V5 %in% SelectedACR_Test,]
bedfile_SelectedACR_control <- ACRs_18cells[ACRs_18cells$V5 %in% SelectedACR_Control_SameNumberwithTest,]

head(bedfile_SelectedACR_Test)
write.table(bedfile_SelectedACR_Test[ , -5], file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs_curated_RandomSelectSNPperACR.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote=F)
write.table(bedfile_SelectedACR_control[ , -5], file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs_curated_RandomSelectSNPperACR.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote=F)
