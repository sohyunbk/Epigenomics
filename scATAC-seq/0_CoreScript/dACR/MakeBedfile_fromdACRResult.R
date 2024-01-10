library("optparse")
library(rlang)


option_list = list(
  make_option(c("--DEGDir"), type="character", 
              help="DEGDir", metavar="character"),
  make_option(c("--DEGFileName"), type="character", 
              help="DEGFileName", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

WDir <- opt$DEGDir
FileEnd <- opt$DEGFileName

#WDir <- "/scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/"
#FileEnd <- ".EdgeRResult_PseudoReplicate_withPromoterRegion.txt"

pattern_string <- paste0(FileEnd, "$")
files <- list.files(path = WDir, pattern = pattern_string, full.names = FALSE)

process_and_write_peaks <- function(subset, file_suffix) {
  peaks <- subset$Peak
  formatted_peaks <- sapply(peaks, function(x) {
    parts <- unlist(strsplit(x, "_"))
    return(paste(parts[1], parts[2], parts[3], sep = "\t"))
  }, USE.NAMES = FALSE)
  writeLines(formatted_peaks, paste0(WDir, name_before_dot, file_suffix))
}

# Write to a file
## logFC >0 --> Bif3 higher
## logFC <0 --> WT higher
for (file in files){
  name_before_dot <- sub("\\..*$", "", file)
  DEGTable <- read.table(file,head=TRUE)
  subset_DEGTable_0.05_Bif3Higher <- DEGTable[DEGTable$FDR < 0.05 & DEGTable$logFC >0 , ]
  subset_DEGTable_0.05_A619Higher <- DEGTable[DEGTable$FDR < 0.05 & DEGTable$logFC <0 , ]

  process_and_write_peaks(subset_DEGTable_0.05_Bif3Higher, ".FDR0.05Bif3Higher.bed")
  process_and_write_peaks(subset_DEGTable_0.05_A619Higher, ".FDR0.05A619Higher.bed")
  
  subset_DEGTable_0.01_Bif3Higher <- DEGTable[DEGTable$FDR < 0.01 & DEGTable$logFC >0 , ]
  subset_DEGTable_0.01_A619Higher <- DEGTable[DEGTable$FDR < 0.01 & DEGTable$logFC <0 , ]
  
  process_and_write_peaks(subset_DEGTable_0.01_Bif3Higher, ".FDR0.01Bif3Higher.bed")
  process_and_write_peaks(subset_DEGTable_0.01_A619Higher, ".FDR0.01A619Higher.bed")
  }



