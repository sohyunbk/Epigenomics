##
A619Meta <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt")
Bif3Meta <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV4_metadata.txt")
head(A619Meta)
AllTn5 <- c(A619Meta$nSites,Bif3Meta$nSites)
length(AllTn5)
mean(AllTn5)
