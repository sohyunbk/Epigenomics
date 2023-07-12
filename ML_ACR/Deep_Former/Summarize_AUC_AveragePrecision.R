RedandantACR <- read.table("/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome_withBigN/test_performance.txt",
                           header=TRUE,row.names="class")
head(RedandantACR)
WithNegative <- read.table("/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative/test_performance.txt",
                           header=TRUE,row.names="class")
head(WithNegative)

WithNegative_SameACRNumber <- read.table("/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_SameNumberNegative/test_performance.txt",
                           header=TRUE,row.names="class")
head(WithNegative_SameACRNumber)

RemoveRedandant <- read.table("/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/DanQ_RemoveRedundantACR/test_performance.txt",
                           header=TRUE,row.names="class")
head(RemoveRedandant)

CellTypeRestrictACR <- read.table("/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/DanQ_CellTypeRestrict/test_performance.txt",
                              header=TRUE,row.names="class")
head(CellTypeRestrictACR)

CombinedTable <- cbind(RedandantACR,WithNegative,WithNegative_SameACRNumber,
                       RemoveRedandant,CellTypeRestrictACR)
head(CombinedTable)
write.table(CombinedTable, "/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/SummaryTable.txt", 
            quote=F, row.names=T, col.names=T, sep="\t")