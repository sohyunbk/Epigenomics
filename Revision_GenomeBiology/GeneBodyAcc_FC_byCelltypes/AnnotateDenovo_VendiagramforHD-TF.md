Getting FC (Bif3/WT) of gene body chromatin acc
================
Sohyun Bang
01 June, 2025

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

### Load Gene Info data

``` r
GeneInfo <- read.delim("../maizev5_data/Zm00001eb.1.fulldata_Curated2.txt", stringsAsFactors = FALSE)
```

### 1) Get Gene body acc for all the genes & normalize the value & calculate logFC

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

GeneXCT_MoreThan50Tn5 <- GeneXCT %>%
  filter_all(all_vars(. > 50))
dim(GeneXCT_MoreThan50Tn5)
```

    ## [1] 18113    29

``` r
GeneNames <- GeneXCT[, 1]
gene_counts_df <- GeneXCT[, -1]

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
  left_join(GeneInfo[, c("gene_model", "locus_symbol")], by = "gene_model")
write.table(FCTable_annotated, file = "GeneBodyACC_FC_Bif3WT_byCelltypes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

CentralZone_Increased <- FCTable_annotated %>%
  filter(IM.OC > 0.5)
```

### 2) Load Homemodomain TFs and match from V3 to V5

``` r
TFlist <- read.delim("./Zma_TF_list.txt", stringsAsFactors = FALSE) ## This data is from https://planttfdb.gao-lab.org/index.php?sp=Zma

## Filter TFs with HD
Temp <- TFlist %>%
  filter(Family %in% c("ZF-HD", "HD-ZIP", "WOX", "HB0PHD", "HB-other"))

V3_V5 <- read.delim("B73v3_to_B73v5.tsv", stringsAsFactors = FALSE)

TFs_merged <- merge(TFlist, V3_V5, by.x = "Gene_ID", by.y = "V3", all.x = TRUE)
TFs_expanded <- TFs_merged %>% ## This is because to make V5 gene id as the key
  separate_rows(V5, sep = ",") %>%
  rename(gene_model=V5) 

TFs_V5Key <- TFs_expanded %>%
  group_by(gene_model) %>%
  summarise(Family = paste(sort(unique(Family)), collapse = ","), .groups = "drop") %>%
  filter(gene_model != "") %>%
  separate_rows(Family, sep = ",") # This is because some gene model has two families

HD_TFs <- TFs_V5Key %>%
    filter(Family %in% c("ZF-HD", "HD-ZIP", "WOX", "HB0PHD", "HB-other")) 

unique_HD_TFs<- unique(HD_TFs$gene_model)
length(unique_HD_TFs)
```

    ## [1] 152

### 3) Load De novo marker gene list

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

### 4) Venn diagram

``` r
common_elements <- intersect(intersect(unique_HD_TFs, CentralZone_Increased$gene_model), Sig_DeNovo_annotated$gene_model)
print(common_elements)
```

    ## [1] "Zm00001eb067310"

``` r
venn_list <- list(HD_TF = unique_HD_TFs, 
                  IncreasedinBif3 = CentralZone_Increased$gene_model,
                  DeNovoCentralzone = Sig_DeNovo_annotated$gene_model)
ggVennDiagram(venn_list, label="count")
```

![](output/fig-Venn%20diagram%20-1.png)<!-- -->
