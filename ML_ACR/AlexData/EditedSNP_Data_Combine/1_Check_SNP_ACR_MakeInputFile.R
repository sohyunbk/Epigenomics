library(dplyr)
library(GenomicRanges)

###Load 18cell types data for bed file later.
CT18 <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/NonRedundantACRs_18Cells.500bp_Sorted.bed")
CT18$V2 <- CT18$V2
CT18$V3 <- CT18$V3
head(CT18)
CT18$Pos <- paste(CT18$V1,CT18$V2,CT18$V3,sep="_")
CT18_combined <- CT18 %>%
  group_by(Pos) %>%
  summarise(V4_combined = paste(V4, collapse = ";")) %>%
  ungroup()
# If you need to merge this with the original data frame to retain all columns:
CT18_final <- CT18 %>%
  distinct(Pos, .keep_all = TRUE) %>%
  left_join(CT18_combined, by = "Pos")
head(CT18_final)

file_path <- "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/NonRedundantACRs_18Cells.500bp_distinctfeatures.txt"
# Read the file as a vector
feature_vector <- readLines(file_path)
# Print the vector to check
print(feature_vector)
feature_map <- setNames(seq_along(feature_vector), feature_vector)

# Function to replace feature names with their corresponding numbers
replace_with_numbers <- function(combined_features) {
  # Split combined features into individual features
  features <- strsplit(combined_features, ";")[[1]]
  # Replace each feature with its corresponding number
  numbers <- sapply(features, function(x) feature_map[x])
  # Combine numbers back into a single string
  paste(numbers, collapse = ";")
}

# Apply the function to each row in CT18_combined$V4_combined
CT18_combined$V4_numbered <- sapply(CT18_combined$V4_combined, replace_with_numbers)

# View the result
head(CT18_combined)


#### ******* Test dataset
SNPFile_test <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs.v2.txt", header=TRUE)
unique_SNPFile_test <- unique(SNPFile_test)
print(nrow(SNPFile_test))
print(nrow(unique_SNPFile_test))
head(SNPFile_test)
dim(SNPFile_test)
element_counts <- table(unique_SNPFile_test$acrID)
head(element_counts)
SNPFile_test[SNPFile_test$snpID == "chr9_99008396",]
head(unique(SNPFile_test$acrID))
length(unique(SNPFile_test$acrID))
## 1) Check the overlap as it's extended with 1000bp
split_pos <- strsplit(unique(SNPFile_test$acrID), "_")
chr <- sapply(split_pos, function(x) x[1])
start <- as.integer(sapply(split_pos, function(x) x[2])) -250
end <- as.integer(sapply(split_pos, function(x) x[3])) +250
gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

overlaps <- findOverlaps(gr)
summary(overlaps)
overlap_pairs <- as.data.frame(overlaps)
overlap_counts <- countOverlaps(gr)
NonOverlappedACRs <- unique(SNPFile_test$acrID)[overlap_counts==1]
head(NonOverlappedACRs)
length(NonOverlappedACRs)
CT18_combined
CT18_combined <- unique(CT18_combined)
NonOverlappedACRs_CTInfo <- CT18_combined[CT18_combined$Pos %in% NonOverlappedACRs, ]

MockData <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR8CT_DanQ/test_data.bed")
length(MockData$V5[0:2401])
NonOverlappedACRs_split <- strsplit(NonOverlappedACRs, "_")
chr <- sapply(NonOverlappedACRs_split, function(x) x[1])
start <- as.integer(sapply(NonOverlappedACRs_split, function(x) x[2])) -250
end <- as.integer(sapply(NonOverlappedACRs_split, function(x) x[3])) +250
BedFile <- data.frame(chr,start,end)
head(BedFile)

BedFile$Strand <- "+"
dim(BedFile)
BedFile$CellTypeTemp <- MockData$V5[0:2401]
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
NonOverlappedACRs <- unique(SNPFile$acrID)[overlap_counts==1]

NonOverlappedACRs_split <- strsplit(NonOverlappedACRs, "_")
chr <- sapply(NonOverlappedACRs_split, function(x) x[1])
start <- as.integer(sapply(NonOverlappedACRs_split, function(x) x[2])) -250
end <- as.integer(sapply(NonOverlappedACRs_split, function(x) x[3])) +250
BedFile <- data.frame(chr,start,end)
dim(BedFile)
BedFile$Strand <- "+"
BedFile$CellTypeTemp <- MockData$V5[0:13414]
dim(BedFile)
write.table(BedFile, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs.v2.curated1000bp.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=F)

head(NonOverlappedACRs)
SNPFile_Filtered <- SNPFile[SNPFile$acrID %in% NonOverlappedACRs,]
dim(SNPFile_Filtered)

UniqueACR <- SNPFile_Filtered %>%
  group_by(acrID) %>%
  sample_n(size = 1) %>%
  ungroup()
write.table(UniqueACR, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs.v2.curated.txt", sep = "\t", row.names = FALSE, quote=F)
