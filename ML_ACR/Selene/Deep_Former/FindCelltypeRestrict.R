## Python command
conda activate r_env
python
import pybedtools
from pybedtools import BedTool

AllPeaks = BedTool("/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks.bed")
AllPeaks.merge(c=4, o="collapse").saveas('/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks_Collapsed.bed')

## R command
library(stringr)
library(ggplot2)

File <- read.table("/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks_Collapsed.bed")
head(File[,c(1:3)])
File$V4[1]
NumberofCellPerACR <- str_count(File$V4, ",")
NumberofCellPerACR[1]
NumberofCellPerACR <- NumberofCellPerACR+1
NumberofCellPerACR[2]

ggplot() + aes(NumberofCellPerACR)+ geom_histogram()
ggsave("/scratch/sb14489/8.ML_ACR/1.InputBed/FrequencyPlot.pdf", width=13, height=10)	


