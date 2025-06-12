A619_Re3<- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/A619_Re3/A619_Re3_loadData.rds")
A619_Re3<- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/A619_Re3/A619_Re3_Tn5Cut1000_Binsize500.rds")
Yinxin <- read.table("//scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/A619/A619_Tn5Cut1000_Binsize500_MinT0.05_MaxT0.05_PC100_FeaturesN4476_k50_res0.9.AfterHarmony.metadata.txt")

Yinxin_Re3 <- Yinxin[Yinxin$sampleID =="A619_Re3",]
dim(Yinxin_Re3)
head(Yinxin_Re3)
median(Yinxin_Re3$total)

Yinxin_Re4 <- Yinxin[Yinxin$sampleID =="A619_Re4",]
dim(Yinxin_Re4)
median(Yinxin_Re4$total)
