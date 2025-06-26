Getting FC (Bif3/WT) of gene body chromatin acc
================
Sohyun Bang
19 June, 2025

    ## Loading required package: limma

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## 
    ## Attaching package: 'ggVennDiagram'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     unite

    ## Loading required package: grid

    ## Loading required package: futile.logger

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ lubridate 1.9.4     ✔ readr     2.1.5
    ## ✔ purrr     1.0.4     ✔ stringr   1.5.1
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()        masks stats::filter()
    ## ✖ dplyr::lag()           masks stats::lag()
    ## ✖ ggVennDiagram::unite() masks tidyr::unite()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

### Load Gene Info data

``` r
GeneInfo <- read.delim("../maizev5_data/Zm00001eb.1.fulldata_Curated2.txt", stringsAsFactors = FALSE)
Gene <- read.delim("../maizev5_data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed",
                   stringsAsFactors = FALSE, header = FALSE)
colnames(Gene)[1:10] <- c("chr", "start", "end", "gene_id", "dot", "strand", "source", "feature", "dot2", "attributes")
Gene$length_kb <- (Gene$end - Gene$start) / 1000
gene_length_vec <- setNames(Gene$length_kb, Gene$gene_id)
```

### 1) Load Homemodomain TFs and match from V3 to V5

``` r
TFlist <- read.delim("./Zma_TF_list.txt", stringsAsFactors = FALSE) ## This data is from https://planttfdb.gao-lab.org/index.php?sp=Zma

## Filter TFs with HD

V3_V5 <- read.delim("B73v3_to_B73v5.tsv", stringsAsFactors = FALSE)

TFs_merged <- merge(TFlist, V3_V5, by.x = "Gene_ID", by.y = "V3", all.x = TRUE)
TFs_expanded <- TFs_merged %>% ## This is because to make V5 gene id as the key
  separate_rows(V5, sep = ",") %>%
  rename(gene_model=V5) 

TFs_V5Key <- TFs_expanded %>%
  group_by(gene_model) %>%
  summarise(Family = paste(sort(unique(Family)), collapse = ","), .groups = "drop") %>%
  filter(gene_model != "") 

TFs_V5Key_re <- TFs_V5Key %>%
  separate_rows(Family, sep = ",") 


HD_TFs <- TFs_V5Key %>%
  separate_rows(Family, sep = ",") %>% # This is because some gene model has two families 
    filter(Family %in% c("ZF-HD", "HD-ZIP", "WOX", "HB-PHD", "HB-other","TALE")) 

#TFs_V5Key %>%
#  filter(gene_model == "Zm00001eb001720")

#HD_TFs %>%
#  filter(gene_model=="Zm00001eb395430")
unique_HD_TFs<- unique(HD_TFs$gene_model)
length(unique_HD_TFs)
```

    ## [1] 184

### 2) Get Gene body acc for all the genes & normalize the value & calculate logFC

``` r
### 1) Get Gene body acc for all the genes
## load gene*cell table.
WT_GeneXCT <- read.table("./A619_AnnV4.GeneBodyACC.byGeneXCT.txt",header=TRUE)
Bif3_GeneXCT <- read.table("./Bif3_AnnV4.GeneBodyACC.byGeneXCT.txt",header=TRUE)
Celltypes <- colnames(WT_GeneXCT)[-1]
colnames(Bif3_GeneXCT)[-1] <- paste("Bif3&", colnames(Bif3_GeneXCT)[-1], sep="")
colnames(WT_GeneXCT)[-1] <- paste("A619&", colnames(WT_GeneXCT)[-1], sep="")

GeneXCT <- merge(WT_GeneXCT, Bif3_GeneXCT, by = "gene")
head(GeneXCT)
```

    ##              gene A619&FloralMeristem_SuppressedBract A619&G2_M A619&IM.OC
    ## 1 Zm00001eb000010                                 178       188        232
    ## 2 Zm00001eb000020                                 463       715        547
    ## 3 Zm00001eb000050                                  18        30         31
    ## 4 Zm00001eb000060                                 210       215        256
    ## 5 Zm00001eb000070                                  53        68         85
    ## 6 Zm00001eb000080                                 535       484        507
    ##   A619&L1 A619&L1atFloralMeristem A619&PhloemPrecursor
    ## 1     168                     414                  448
    ## 2     402                     930                 1649
    ## 3      22                      58                   48
    ## 4     185                     364                  463
    ## 5      47                     137                  155
    ## 6     479                     927                 1175
    ##   A619&ProcambialMeristem_ProtoXylem_MetaXylem
    ## 1                                          207
    ## 2                                          771
    ## 3                                           34
    ## 4                                          254
    ## 5                                           61
    ## 6                                          620
    ##   A619&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
    ## 1                                                         71
    ## 2                                                        238
    ## 3                                                         29
    ## 4                                                         93
    ## 5                                                         40
    ## 6                                                        241
    ##   A619&SPM.base_SM.base A619&Unknown_lowFRiP A619&Unknown_Sclerenchyma
    ## 1                   179                  347                       287
    ## 2                   402                  669                       577
    ## 3                    36                   66                        26
    ## 4                   180                  251                       267
    ## 5                    45                  167                       100
    ## 6                   462                  593                       720
    ##   A619&Unknown1 A619&Unknown2 A619&XylemParenchyma_PithParenchyma
    ## 1           272           248                                 263
    ## 2           589           435                                 906
    ## 3            27            29                                  46
    ## 4           228           181                                 355
    ## 5            87            93                                  90
    ## 6           740           565                                 913
    ##   Bif3&FloralMeristem_SuppressedBract Bif3&G2_M Bif3&IM.OC Bif3&L1
    ## 1                                 160        89         40      95
    ## 2                                 322       307        158     297
    ## 3                                  21        11          6      12
    ## 4                                 136        86         52     115
    ## 5                                  37        25         39      56
    ## 6                                 394       271        146     268
    ##   Bif3&L1atFloralMeristem Bif3&PhloemPrecursor
    ## 1                     148                  260
    ## 2                     455                  946
    ## 3                      34                   30
    ## 4                     147                  237
    ## 5                      85                   87
    ## 6                     411                  621
    ##   Bif3&ProcambialMeristem_ProtoXylem_MetaXylem
    ## 1                                          133
    ## 2                                          379
    ## 3                                           17
    ## 4                                          131
    ## 5                                           44
    ## 6                                          281
    ##   Bif3&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma
    ## 1                                                         50
    ## 2                                                        171
    ## 3                                                          8
    ## 4                                                         55
    ## 5                                                         22
    ## 6                                                        173
    ##   Bif3&SPM.base_SM.base Bif3&Unknown_lowFRiP Bif3&Unknown_Sclerenchyma
    ## 1                    74                  137                       159
    ## 2                   210                  334                       321
    ## 3                     4                   23                        11
    ## 4                    70                  106                       156
    ## 5                    24                   62                        66
    ## 6                   255                  272                       425
    ##   Bif3&Unknown1 Bif3&Unknown2 Bif3&XylemParenchyma_PithParenchyma
    ## 1           148           140                                 129
    ## 2           349           294                                 401
    ## 3            15             3                                  27
    ## 4           135           102                                 146
    ## 5            39            41                                  48
    ## 6           377           314                                 356

``` r
## normalization
GeneNames <- GeneXCT[, 1]
gene_counts_df <- GeneXCT[, -1]

## Set up the cut off to 50!
#head(GeneXCT)
#dim(GeneXCT)


GeneXCT_long <- GeneXCT %>%
  pivot_longer(cols = -gene, names_to = "CellType", values_to = "Tn5")

ggplot(GeneXCT_long, aes(x = Tn5)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of Tn5 Values",
       x = "Tn5 Signal",
       y = "Density") +
  theme_minimal() 
```

![](output/fig-CheckTn5CountsDistributions%7D-1.png)<!-- -->

``` r
ggplot(GeneXCT_long, aes(x = Tn5)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of Tn5 Values",
       x = "Tn5 Signal",
       y = "Density") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 500)) 
```

![](output/fig-CheckTn5CountsDistributions%7D-2.png)<!-- -->

``` r
#GeneXCT %>%
#  filter(gene == "Zm00001eb067310")
#176
#GeneXCT %>%
#  filter(gene == "Zm*00001eb395430")
```

Normalization - TPM considering gene length - due to the cut off.

``` r
## Here only filter the genes that are less than 50 counts in central zone..
GeneXCT_MoreThan50Tn5 <- subset(GeneXCT, `A619&IM.OC` > 50 & `Bif3&IM.OC` > 50)
expr <- GeneXCT_MoreThan50Tn5
matched_genes <- intersect(expr$gene, names(gene_length_vec))
expr <- expr[expr$gene %in% matched_genes, ]
gene_lengths_kb <- gene_length_vec[expr$gene]
stopifnot(all(expr$gene == names(gene_lengths_kb))) 
count_matrix <- as.matrix(expr[, -1]) 

calculateTPM <- function(counts, gene_lengths_kb) {
  rate <- counts / gene_lengths_kb
  tpm <- sweep(rate, 2, colSums(rate), FUN = "/") * 1e6
  return(tpm)
}

tpm_matrix <- calculateTPM(count_matrix, gene_lengths_kb)
rownames(tpm_matrix) <- expr$gene

#check ZmWUS1
tpm_matrix["Zm00001eb067310", ] # itss tpm value is 90.60476
```

    ##                        A619&FloralMeristem_SuppressedBract 
    ##                                                   52.49291 
    ##                                                  A619&G2_M 
    ##                                                   35.33419 
    ##                                                 A619&IM.OC 
    ##                                                   83.17461 
    ##                                                    A619&L1 
    ##                                                   29.53491 
    ##                                    A619&L1atFloralMeristem 
    ##                                                   47.71927 
    ##                                       A619&PhloemPrecursor 
    ##                                                   48.76834 
    ##               A619&ProcambialMeristem_ProtoXylem_MetaXylem 
    ##                                                   39.99035 
    ## A619&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma 
    ##                                                   37.91976 
    ##                                      A619&SPM.base_SM.base 
    ##                                                   67.32539 
    ##                                       A619&Unknown_lowFRiP 
    ##                                                   41.30202 
    ##                                  A619&Unknown_Sclerenchyma 
    ##                                                   38.35709 
    ##                                              A619&Unknown1 
    ##                                                   45.01292 
    ##                                              A619&Unknown2 
    ##                                                   32.69206 
    ##                        A619&XylemParenchyma_PithParenchyma 
    ##                                                   30.90417 
    ##                        Bif3&FloralMeristem_SuppressedBract 
    ##                                                   99.61386 
    ##                                                  Bif3&G2_M 
    ##                                                   86.19118 
    ##                                                 Bif3&IM.OC 
    ##                                                  116.34624 
    ##                                                    Bif3&L1 
    ##                                                   65.53095 
    ##                                    Bif3&L1atFloralMeristem 
    ##                                                  104.00986 
    ##                                       Bif3&PhloemPrecursor 
    ##                                                   82.16736 
    ##               Bif3&ProcambialMeristem_ProtoXylem_MetaXylem 
    ##                                                   60.60903 
    ## Bif3&ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma 
    ##                                                   73.65183 
    ##                                      Bif3&SPM.base_SM.base 
    ##                                                  119.43668 
    ##                                       Bif3&Unknown_lowFRiP 
    ##                                                   90.60476 
    ##                                  Bif3&Unknown_Sclerenchyma 
    ##                                                   87.45195 
    ##                                              Bif3&Unknown1 
    ##                                                  103.05610 
    ##                                              Bif3&Unknown2 
    ##                                                   70.32274 
    ##                        Bif3&XylemParenchyma_PithParenchyma 
    ##                                                   67.77750

``` r
Genes_higherTPMthanBif3CZ <- rownames(tpm_matrix[tpm_matrix[,"Bif3&IM.OC"] >= 90.50476, ])
 

GeneNames <- GeneXCT_MoreThan50Tn5[, 1]
gene_counts_df <- GeneXCT_MoreThan50Tn5[, -1]

cpm_data <- cpm(DGEList(counts = gene_counts_df), log = FALSE)
QNorm <- normalize.quantiles(as.matrix(gene_counts_df))
rownames(QNorm) <- GeneNames
colnames(QNorm) <- colnames(GeneXCT[,-1])
#head(QNorm)
#dim(QNorm)

# Build fold change table
FCTable <- as.data.frame(sapply(Celltypes, function(Celltypes) {
  bif3_col <- paste0("Bif3&", Celltypes)
  A619_col <- paste0("A619&", Celltypes)
  log2(QNorm[, bif3_col] / QNorm[, A619_col])
}))

#Matching with gene names
FCTable_with_id <- FCTable %>%
  tibble::rownames_to_column(var = "gene_model")

FCTable_annotated <- FCTable_with_id %>%
  left_join(GeneInfo[, c("gene_model", "locus_symbol")], by = "gene_model") %>%
  left_join(TFs_V5Key[, c("gene_model", "Family")], by = "gene_model")

write.table(FCTable_annotated, file = "GeneBodyACC_FC_Bif3WT_byCelltypes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

CentralZone_Increased <- FCTable_annotated %>%
  filter(IM.OC > 0.5)

CentralZone_Increased_morethanBif3CZ_HD <- CentralZone_Increased %>%
  filter(gene_model %in% Genes_higherTPMthanBif3CZ) %>%
  filter(gene_model %in% unique_HD_TFs)

write.table(CentralZone_Increased_morethanBif3CZ_HD, file = "CentralZone_Increased_morethanBif3CZ_HD.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#CentralZone_Increased_morethanBif3CZ %>%
#  filter(gene_model == "Zm00001eb067310")
```

### 3) Load De novo marker gene list \<\<- its not used

``` r
DeNovo <- read.delim("A619_v4.IM.OC_deseq_2_results.tsv", stringsAsFactors = FALSE)
Sig_DeNovo <- DeNovo[!is.na(DeNovo$padj) & 
                     !is.na(DeNovo$log2FoldChange) & 
                     DeNovo$pvalue < 0.01 & 
                     DeNovo$log2FoldChange > 0, ]

Sig_DeNovo_annotated <- Sig_DeNovo %>%
  rename(gene_model = gene_name) %>%
  left_join(GeneInfo[, c("gene_model", "locus_symbol")], by = "gene_model")
write.table(Sig_DeNovo_annotated, file = "Sig_DeNovo_p0.01_genesymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

### 4) barplot

``` r
#common_elements <- intersect(intersect(unique_HD_TFs, CentralZone_Increased$gene_model), Sig_DeNovo_annotated$gene_model)
common_elements <- intersect(unique_HD_TFs, CentralZone_Increased$gene_model)

temp <- TFs_V5Key_re[TFs_V5Key_re$gene_model %in% CentralZone_Increased$gene_model,]
full_table <- table(TFs_V5Key_re$Family)
temp_table <- table(temp$Family)

all_families <- names(full_table)
temp_vector <- setNames(rep(0, length(all_families)), all_families)
temp_vector[names(temp_table)] <- temp_table

ratio <- temp_vector / full_table * 100

ratio_df <- data.frame(
  Family = names(ratio),
  Temp_Count = temp_vector,
  Full_Count = as.numeric(full_table),
  Ratio = as.numeric(ratio)
)
#print(ratio_df)
ratio_df_filtered <- ratio_df %>% 
  filter(Ratio > 0)

# Plot
ggplot(ratio_df_filtered, aes(x = fct_reorder(Family, Ratio), y = Ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Temp_Count), 
            hjust = -0.2, size = 4) +  # Adjust size and position as needed
  coord_flip() +
  labs(title = "Ratio of Temp_Count to Full_Count by TF Family (Non-zero Only)",
       x = "Transcription Factor Family",
       y = "Ratio (%)") +
  theme_minimal(base_size = 14) +
  ylim(0, max(ratio_df_filtered$Ratio) * 1.1)  
```

![](output/fig-barplot%20-1.png)<!-- --> \### 5) Venn diagram

``` r
venn_list <- list(HD_TF = unique_HD_TFs, 
                  IncreasedinBif3 = CentralZone_Increased$gene_model)
#ggVennDiagram(venn_list, label="count")
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("HD_TF", "IncreasedinBif3"),
  filename = NULL,  # Plot to R graphics window
  output = TRUE,
  imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)
grid::grid.draw(venn.plot)
```

![](output/fig-Venn%20diagram%20-1.png)<!-- -->
