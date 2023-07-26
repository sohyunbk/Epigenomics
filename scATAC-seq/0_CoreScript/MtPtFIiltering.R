### MtPt Density plot

Run <-function(SampleName) {

obj <- readRDS(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/",
                      SampleName, "/", SampleName,"_Tn5Cut1000_Binsize500.rds"))
str(obj)

head(obj$meta)

library(tidyverse)
library(ggplot2)
library(grid)

Cutoffcell<- sum(obj$meta$pPtMt < 0.05)

# Bin size control + color palette
# Create the ggplot with the desired aesthetics
p <- ggplot(obj$meta, aes(x = log10nSites, y = pPtMt)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  xlab("Tn5 integration sites per barcode (log10)") +
  ylab("Organelle Ratio") +
  theme_bw() +
  theme(legend.text = element_blank())

# Add the line across y = 0.25
p <- p + geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")

# Calculate the position to place the custom text on the right top corner
x_pos <- max(obj$meta$log10nSites) - 0.1
y_pos <- max(obj$meta$pPtMt) - 0.03

# Add the custom text on the right top corner
p <- p + annotate(
  "text",
  x = x_pos, y = y_pos,
  label = paste0("# Cell = ", Cutoffcell),
  hjust = 1, vjust = 1,
  size = 4, color = "black"
)

# Save the plot as a PDF
ggsave(
  paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/",
   SampleName,"/",SampleName,"Org.pdf"),
  p,
  width = 7, height = 5,
  units = "in", dpi = 300
)
}

Run("A619_Re3")
Run("A619_Re4")

Run("bif3_Re3")
Run("bif3_Re4")


Org_Frip <-function(SampleName) {
  
  obj <- readRDS(paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/",
                        SampleName, "/", SampleName,"_Tn5Cut1000_Binsize500.rds"))
  str(obj)
  
  head(obj$meta)
  
  library(tidyverse)
  library(ggplot2)
  library(grid)
  
  Cutoffcell<- sum(obj$meta$pPtMt < 0.05)
  Meta <- obj$meta
  Meta$FRiP <- Meta$acr/Meta$total
  head(Meta)
  # Bin size control + color palette
  ggplot(Meta, aes(x = FRiP, y = pPtMt)) +
    geom_bin2d(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    xlab("FRiP") +
    ylab("Organelle Ratio") +
    theme_bw() +
    theme(legend.text = element_blank())
  
  ggsave(
    paste0("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/",
           SampleName, "/", SampleName,
           "_org_FRiP.pdf"),width = 7, height = 5)
    
}

Org_Frip("A619_Re3")
Org_Frip("A619_Re4")

Org_Frip("bif3_Re3")
Org_Frip("bif3_Re4")