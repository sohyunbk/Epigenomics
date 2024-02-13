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
dACR <- (ACR_ClosestGene[ACR_ClosestGene$FDR < 0.05,])
table(dACR$TF)
dim(dACR)
table(dACR$TFFamily)

dim(dACR)

write.table(dACR, 
            "/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.XSTREME/AnnV4/IM-OC.FDR0.05Bif3Higher.ControlfromIntergenicAllSameCTPeaks.XSTREME/dACR_withTAATInfo.txt", 
            quote=F, row.names=F, col.names=T, sep="\t")

