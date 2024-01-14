library(dplyr)

ACRs_18cells <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/Seedling_18Celltypes.500.RestrictACR18CT.bed",sep="\t", header=FALSE, stringsAsFactors=FALSE)
ACRs_pos <- paste(ACRs_18cells$V1, ACRs_18cells$V2, ACRs_18cells$V3, sep="_")
ACRs_pos_NonRe <- unique(ACRs_pos)
length(ACRs_pos_NonRe)
ACRs <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/AllACRs_130539.500bp.Sorted.bed",sep="\t", header=FALSE, stringsAsFactors=FALSE)
ACRs_pos1 <- paste(ACRs$V1, ACRs$V2, ACRs$V3, sep="_")
ACRs_pos1NonRe <- unique(ACRs_pos1)
length(ACRs_pos1NonRe)
head(ACRs_pos1NonRe)

#### 1) Remove redundant lines in snp files
SNPFile_control <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs.txt", header=TRUE)
unique_SNPFile_control <- unique(SNPFile_control)
print(nrow(SNPFile_control))
print(nrow(unique_SNPFile_control))
element_counts <- table(unique_SNPFile_control$snpID)
head(element_counts)
duplicates <- names(element_counts[element_counts > 1])
head(unique_SNPFile_control)
length(unique(unique_SNPFile_control$acrID))
write.table(unique_SNPFile_control, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs_curated.txt", sep = "\t", row.names = FALSE, quote=F)
set.seed(123) # Setting a seed for reproducibility of random selection
randomly_selected_control <- unique_SNPFile_control %>%
  group_by(acrID) %>%
  sample_n(size = 1)
dim(randomly_selected_control)
head(randomly_selected_control)
unique_SNPFile_control[c(1:20),]
randomly_selected_control[randomly_selected_control$acrID =="chr1_167642_168142",]
print(randomly_selected_control[randomly_selected_control$acrID =="chr1_187672_188172",]$snpID)
write.table(randomly_selected_control, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs_curated_RandomSelectSNPperACR.txt", sep = "\t", row.names = FALSE, quote=F)

#split_snp <- strsplit(unique_SNPFile_control$snpID, "_")
# Extract the parts you want
#snp_data <- data.frame(
#  chr = sapply(split_snp, "[[", 1),
#  pos = sapply(split_snp, "[[", 2),
#  alt = unique_SNPFile_control$alt
#)
# Print the first few rows of the new table
#head(snp_data)
#write.table(snp_data, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs_curated_pythonInput.txt", sep = "\t", row.names = FALSE,col.names=FALSE, quote=F)



SNPFile_test <- read.table("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs.txt", header=TRUE)
unique_SNPFile_test <- unique(SNPFile_test)
print(nrow(SNPFile_test))
print(nrow(unique_SNPFile_test))
element_counts <- table(unique_SNPFile_test$snpID)
head(element_counts)
duplicates <- names(element_counts[element_counts > 1])
SNPFile_test[SNPFile_test$snpID == "chr9_99008396",]

merged_SNPFile_test <- aggregate(. ~ acrID + snpID + ref + alt, SNPFile_test, function(x) if (is.numeric(x)) max(x) else x)
merged_SNPFile_test[merged_SNPFile_test$snpID == "chr9_99008396",]
nrow(merged_SNPFile_test)
head(merged_SNPFile_test)
length(unique(merged_SNPFile_test$acrID))
tail(merged_SNPFile_test)

write.table(merged_SNPFile_test, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs_curated.txt", sep = "\t", row.names = FALSE, quote=F)

randomly_selected_test <- merged_SNPFile_test %>%
  group_by(acrID) %>%
  sample_n(size = 1)

write.table(randomly_selected_test, file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs_curated_RandomSelectSNPperACR.txt", sep = "\t", row.names = FALSE, quote=F)

#### 2) Make Input bed file! :)
### Check common ACRs
head(randomly_selected_control$acrID)
head(ACRs_pos1NonRe) 
head(randomly_selected_control$acrID)
length(intersect(randomly_selected_control$acrID, ACRs_pos_NonRe))
length(intersect(merged_SNPFile_test$acrID, ACRs_pos_NonRe))
SelectedACR_Control <- intersect(randomly_selected_control$acrID, ACRs_pos_NonRe)
SelectedACR_Test <- intersect(merged_SNPFile_test$acrID, ACRs_pos_NonRe)
SelectedACR_Control_SameNumberwithTest <- sample(SelectedACR_Control, size = length(intersect(merged_SNPFile_test$acrID, ACRs_pos_NonRe)))
head(SelectedACR_Control_SameNumberwithTest)

head(ACRs_18cells)
ACRs_18cells$V5 <- paste(ACRs_18cells$V1, ACRs_18cells$V2, ACRs_18cells$V3, sep="_")
bedfile_SelectedACR_Test <- ACRs_18cells[ACRs_18cells$V5 %in% SelectedACR_Test,]
bedfile_SelectedACR_control <- ACRs_18cells[ACRs_18cells$V5 %in% SelectedACR_Control_SameNumberwithTest,]

head(bedfile_SelectedACR_Test)
write.table(bedfile_SelectedACR_Test[ , -5], file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs_curated_RandomSelectSNPperACR.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote=F)
write.table(bedfile_SelectedACR_control[ , -5], file = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs_curated_RandomSelectSNPperACR.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote=F)
