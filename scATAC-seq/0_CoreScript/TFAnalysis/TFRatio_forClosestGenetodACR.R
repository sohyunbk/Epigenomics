library(GenomicRanges)

## Load Fimo result
Fimo_TAAT <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/fimo_out_2/fimo.tsv",header=TRUE)
dim(Fimo_TAAT)
head(Fimo_TAAT)
Fimo_TAAT_gr <- makeGRangesFromDataFrame(Fimo_TAAT, seqnames.field = "sequence_name", start.field = "start", end.field = "stop", strand.field = "strand", keep.extra.columns = TRUE)

## Load Closest dACR file!
ACR_ClosestGene <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/IM-OC.EdgeRResult_PseudoReplicate_withPromoterRegion_NearestGENEINFO.txt"
                               , fill = TRUE,header=TRUE)
head(ACR_ClosestGene)
dim(ACR_ClosestGene)

#dACR_Bif3Higher <- (ACR_ClosestGene[ACR_ClosestGene$logFC > 0 & ACR_ClosestGene$FDR < 0.05,])
split_names <- strsplit(ACR_ClosestGene$Peak, split = "_")
chromosomes <- sapply(split_names, function(x) x[1])
starts <- as.integer(sapply(split_names, function(x) x[2]))
ends <- as.integer(sapply(split_names, function(x) x[3]))
ACR_gr <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))

overlaps <- findOverlaps(ACR_gr, Fimo_TAAT_gr)
query_hits_indices <- queryHits(overlaps)
length(query_hits_indices)
length(unique(query_hits_indices))
ACR_ClosestGene$TAAT <- "None"
ACR_ClosestGene[query_hits_indices,]$TAAT <- "TAAT"
tail(ACR_ClosestGene)
###

TFFile <- read.csv("/scratch/sb14489/3.scATAC/2.Maize_ear/11.ComprehensiveAnalysis/1.TF_fordACR/data.csv")
head(TFFile)
TFFamilyProfile <- table(TFFile$family)

ACR_ClosestGene$TF <- "NoTF"
ACR_ClosestGene$TFFamily <- "NoTF"

ACR_ClosestGene[ACR_ClosestGene$gene_model %in% TFFile$gene.ID,]$TF <- "TF"

# Ensure the gene_model column is character type for accurate matching
ACR_ClosestGene$gene_model <- as.character(ACR_ClosestGene$gene_model)
# Create a named vector of families from TFFile, with gene.IDs as names
family_vector <- setNames(TFFile$family, TFFile$gene.ID)
# Use gene_model to look up corresponding family, default to "NoTF" if not found
ACR_ClosestGene$TFFamily <- family_vector[ACR_ClosestGene$gene_model]
# Replace NAs with "NoTF" in case there are gene_model values without matching TF families
ACR_ClosestGene$TFFamily[is.na(ACR_ClosestGene$TFFamily)] <- "NoTF"
# Now ACR_ClosestGene has a new column "TFFamily" with the matched family or "NoTF"
tail(ACR_ClosestGene)
head(ACR_ClosestGene)

dACR <- (ACR_ClosestGene[ACR_ClosestGene$FDR < 0.05,])
table(dACR$TF)
dim(dACR)
table(dACR$TFFamily)

dim(dACR)

write.table(dACR, 
            "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/", 
            quote=F, row.names=F, col.names=T, sep="\t")

### Random sampling from ACR.
set.seed(123) # Optional: for reproducibility
random_numbers <- sample(1:37318, 3194, replace = FALSE)
Random_ACR <- ACR_ClosestGene[random_numbers,]
table(Random_ACR$TF)
table(Random_ACR$TFFamily)

dACR <- (ACR_ClosestGene[ACR_ClosestGene$FDR < 0.05,])

dACR_Bif3Higher <- (ACR_ClosestGene[ACR_ClosestGene$logFC > 0 & ACR_ClosestGene$FDR < 0.05,])
dACR_Bif3Higher_withTAAT <- (ACR_ClosestGene[ACR_ClosestGene$logFC > 0 & 
                                            ACR_ClosestGene$FDR < 0.05 & ACR_ClosestGene$TAAT =="TAAT",])

dACR_A619Higher <- (ACR_ClosestGene[ACR_ClosestGene$logFC < 0 & ACR_ClosestGene$FDR < 0.05,])

head(dACR_Bif3Higher_withTAAT)
dim(dACR_Bif3Higher)
table(dACR_Bif3Higher$TF)
table(dACR_A619Higher$TF)
table(dACR_Bif3Higher_withTAAT$TF)
table(dACR_Bif3Higher_withTAAT$TFFamily)
TFFamilyProfile[table(dACR_Bif3Higher_withTAAT$TFFamily)]
TFNumber_TAAT <- table(dACR_Bif3Higher_withTAAT$TFFamily)
TFNumber_TAAT <- TFNumber_TAAT[names(TFNumber_TAAT) != "NoTF"]
TFFamilyProfile[names(TFNumber_TAAT)]
TFNumber_TAAT/TFFamilyProfile[names(TFNumber_TAAT)]

TFNumber_dACR <- table(dACR$TFFamily)
TFNumber_dACR <- TFNumber_dACR[names(TFNumber_dACR) != "NoTF"]
TFFamilyProfile[names(TFNumber_dACR)]
TFNumber_dACR/TFFamilyProfile[names(TFNumber_dACR)]

dACR_Bif3Higher$locus_symbol
dACR_Bif3Higher_withTAAT$locus_symbol
dACR_Bif3Higher_withTAAT$Peak

