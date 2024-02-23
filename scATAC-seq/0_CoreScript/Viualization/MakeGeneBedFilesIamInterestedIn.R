## Load bed File
## The file should be: chr	start	end	geneID	name	type
GeneBedLoad <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr_AddGeneSymbol.bed",header=FALSE)
head(GeneBedLoad)
GeneBed <- GeneBedLoad[grep("^(wox|ZmWus|arftf|pin)", GeneBedLoad$V10, ignore.case = TRUE), ]
GeneBed <- rbind(GeneBed,GeneBedLoad[GeneBedLoad$V4=="Zm00001eb433010",])
GeneBed <- rbind(GeneBed,GeneBedLoad[GeneBedLoad$V4=="Zm00001eb067310",])

OutFile <- data.frame(chr =GeneBed$V1,
                      start=GeneBed$V2,
                      end=GeneBed$V3,
                      geneID=GeneBed$V4,
                      name=GeneBed$V10,
                      type=GeneBed$V10)
OutFile$type <- gsub("\\d", "", OutFile$type)
write.table(OutFile,
            "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/ARF_WOX.txt",
            quote=F, row.names=F, col.names=T, sep="\t")
