Bif3Higher <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_FDR0.01_Bif3Higher_ControlIMOC_XSTREME/streme_out/sequences.tsv", header = TRUE, sep = "\t")
head(Bif3Higher)
tail(Bif3Higher)
dim(Bif3Higher)
Bif3Higher$seq_Class == "fp"

dACR_Bif3Higher <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_BiggerPeaks_AllIntergenic_SeedOn/IM-OC_FDR.0.01_Bif3Higher.Bed")
head(dACR_Bif3Higher)
All <- paste0(dACR_Bif3Higher$V1,":",dACR_Bif3Higher$V2,"-",dACR_Bif3Higher$V3)
head(All)
length(All)
STREME_4 <- Bif3Higher$seq_ID[Bif3Higher$motif_ALT_ID =="STREME-4"]
STREME_5 <- Bif3Higher$seq_ID[Bif3Higher$motif_ALT_ID =="STREME-5"]

length(STREME_4)
length(unique(STREME_4))
length(STREME_5)
length(unique(STREME_5))

STREME_4_Filtered <- intersect(STREME_4,All)
length(STREME_4_Filtered)
STREME_5_Filtered <- intersect(STREME_5,All)
length(STREME_5_Filtered)

intersect(STREME_4_Filtered,STREME_5_Filtered)

length(intersect(STREME_4,All))

total_cards <- 786
cards_round1 <- 388
cards_round2 <- 430
overlap_cards <- 226

# Calculate the overall probability of picking an overlapped card
p_overlap <- (cards_round1 / total_cards) * (cards_round2 / total_cards)

# Calculate the probability of observing 226 or more overlapped cards under the null hypothesis
p_value <- pbinom(overlap_cards, size = total_cards, prob = p_overlap, lower.tail = FALSE)

# Print the p-value
p_value


library(ggvenn)
library(ggplot2)
p <- ggvenn(list(STREME_4_ATWAT = STREME_4, STREME_5_CATG = STREME_5))
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/15.MEME_Motif/IM_OC_dACR_FDR0.01_Bif3Higher_ControlIMOC_XSTREME")
ggsave("Bif3Higher_STREME4_STREME5.pdf", plot = p, width = 6, height = 4) 