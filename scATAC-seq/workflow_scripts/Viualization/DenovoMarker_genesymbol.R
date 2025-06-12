IMOCDenovo <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV4/A619/A619_v4.IM-OC.upregulated_genes.deseq2_output.tsv",header=TRUE)
head(IMOCDenovo)
dim(IMOCDenovo)
GeneInfo <- read.table("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.Curated.txt",
                       fill = TRUE,header=TRUE)
GeneID_GeneSymbol <- data.frame(GeneInfo$gene_model,GeneInfo$locus_symbol)
head(IMOCDenovo)
merged_table <- merge(IMOCDenovo, GeneID_GeneSymbol, by.x = "gene_name", by.y = "GeneInfo.gene_model", all.x = TRUE)
head(merged_table)
dim(merged_table)
#merged_table$GeneInfo.locus_symbol
merged_table[grep("^wox", merged_table$GeneInfo.locus_symbol),]
merged_table[grep("^knox", merged_table$GeneInfo.locus_symbol),]
merged_table[grep("^arf", merged_table$GeneInfo.locus_symbol),]
