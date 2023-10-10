### This script is using Bed file. To use the txt file the matrix, it's different code!

library(ggplot2)

setwd("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/")

AllBed <- read.table("AllACRs_130539.500bp.Sorted.bed",header=FALSE)
head(AllBed)
Selected18Cells <- c("abaxial_bundle_sheath","adaxial_leaf_primordia", "bundle_sheath", "cortex", "dividing_leaf_primordia", "ground_meristem", "guard_mother_cell", 
                     "hypodermal_sclerenchyma", "L1_leaf_primordia_boundary", "leaf_primordia", "mesophyll", "mesophyll_precursors", "phloem_SE_procambial_precursors", "pith_parenchyma", 
                     "procambial_meristem", "protodermal_cell", "protophloem_SE", "xylem_parenchyma")
#CellType18_Past <- read.table("Seedling_18Celltypes.500.bed.sorted",header=FALSE)
CellType18 <- AllBed[AllBed$V4 %in% Selected18Cells,]
head(CellType18)
dim(CellType18)
unique(CellType18$celltype)
#dim(CellType18_Past)
write.table(CellType18[,c(1:4)],
            file ="NonRedundantACRs_18Cells.500bp.bed",
            col.names=FALSE, sep="\t",
            quote=F, row.names=F)

DrawHistogram_forCelltypeRestrictACR <- function(CellType18,SavingName){
  colnames(CellType18) <- c("chr","start","end","celltype")
  CellType18$Pos <- paste(CellType18$chr,CellType18$start,CellType18$end,sep="_")
  head(CellType18)
  
  freq_df <- table(CellType18$Pos)
  freq_df <- as.data.frame(freq_df)
  
  head(freq_df)
  ggplot(freq_df, aes(x=Freq)) + 
    geom_histogram(binwidth=1)  +
    labs(title = " ", x = "Cell type number", y = "Count")
  
  ggsave(paste0(SavingName,".pdf"),width=7,height=6)
  return(freq_df)
}
############
freq_df_18Cells <- DrawHistogram_forCelltypeRestrictACR(CellType18,
                                     "Seedling_18Cells_Histogram_forACRs")
freq_df_AllCells <- DrawHistogram_forCelltypeRestrictACR(AllBed,
                                      "Seedling_AllCells_Histogram_forACRs")
##############

Filter_Cells <- function(CellType18,freq_df,CutOff){
  FilteredACRs <- freq_df$Var1[freq_df$Freq <= CutOff]
  FilteredBed <- CellType18[CellType18$Pos %in% FilteredACRs,]
  write.table(FilteredBed[,c(1:4)],
              file =paste0("Seedling_18Celltypes.500.RestrictACR",CutOff,"CT.bed"),
              col.names=FALSE, sep="\t",
              quote=F, row.names=F)
}

for (cutoff in c(2:18)){
  print(cutoff)
  Filter_Cells(CellType18,freq_df_18Cells,cutoff)
}
