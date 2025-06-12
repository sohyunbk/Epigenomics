A619_Meta <- read.table("/Users/sohyun/Documents/2.SingleCellATAC/3.Data/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt")
head(A619_Meta)
table(A619_Meta$sampleID)
A619_Meta_Re1<-A619_Meta[A619_Meta$sampleID == "A619_Re1",]
A619_Meta_Re2<-A619_Meta[A619_Meta$sampleID == "A619_Re2",]

tail(A619_Meta_Re1)
dim(A619_Meta_Re2)
mean(A619_Meta_Re1$total)
mean(A619_Meta_Re2$total)

mean(A619_Meta_Re1$pTSS)
mean(A619_Meta_Re2$pTSS)

mean(A619_Meta_Re1$FRiP)
mean(A619_Meta_Re2$FRiP)

mean(A619_Meta_Re1$pPtMt)
mean(A619_Meta_Re2$pPtMt)


Bif3_Meta <- read.table("/Users/sohyun/Documents/2.SingleCellATAC/3.Data/6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV4_metadata.txt")
table(Bif3_Meta$sampleID)
tail(Bif3_Meta)
Bif3_Meta_Re1<-Bif3_Meta[Bif3_Meta$sampleID == "bi3_Re1",]
Bif3_Meta_Re2<-Bif3_Meta[Bif3_Meta$sampleID == "bi3_Re2",]

mean(Bif3_Meta_Re1$total)
mean(Bif3_Meta_Re2$total)

mean(Bif3_Meta_Re1$pTSS)
mean(Bif3_Meta_Re2$pTSS)

mean(Bif3_Meta_Re1$FRiP)
mean(Bif3_Meta_Re2$FRiP)

mean(Bif3_Meta_Re1$pPtMt)
mean(Bif3_Meta_Re2$pPtMt)

Total <- rbind(A619_Meta,Bif3_Meta)
mean(Total$total)
