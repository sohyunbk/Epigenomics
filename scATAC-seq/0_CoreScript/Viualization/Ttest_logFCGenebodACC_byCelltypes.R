###### 
library(tidyr)
library(dplyr)
library(tibble)


LogFC <- read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/11.dACR_Character/2.dACR_GeneBodyACC/logFCNeartGene_TAATdACR.txt")
head(LogFC)

LogFC_cleaned <- LogFC %>% 
  select(-Unknown_Sclerenchyma, -Unknown_lowFRiP, -G2_M, -Unknown1,-Unknown2, -XylemParenchyma_PithParenchyma )

LogFC_long <- LogFC_cleaned %>%
  rownames_to_column("Gene") %>%
  gather(key = "CellType", value = "logFC", -Gene)

pairwise_results <- pairwise.t.test(LogFC_long$logFC, LogFC_long$CellType, p.adjust.method = "BH")
print(pairwise_results)
