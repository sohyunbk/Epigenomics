## This script is to draw the QC plot.
## This is because the QC script for original samples do not have
library(ggplot2)
library(viridis)
library("optparse")

## Should it be before QC or after QC? Let's do before QC.-from all bed file.

option_list = list(
  make_option(c("--SampleName"), type="character",
              help="SampleName", metavar="character"),
  make_option(c("--Path"), type="character",
              help="Path")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Samples <- c("1_A619",
#             "1_A619_2",
#             "2_rel2",
#             "2_rel2_2",
#             "3_bif3",
#             "3_bif3_2")

#Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/"

setwd(paste0(opt$Path,opt$SampleName))
print(sample)
obj <- readRDS(paste0(sample,"_Tn5Cut1000_Binsize500.rds"))
Cutoffcell <- sum(obj$meta$pPtMt < 0.05)[TRUE]

colors <- c(00,00, rev(viridis(200, option = "mako"))[2:150])
values <- c(0, seq(0.01, 1, length.out = 200))
p <- ggplot(obj$meta, aes(x = log10nSites, y = pPtMt)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 200) +
  scale_fill_gradientn(colors = colors, values = values) +  # Use the custom color scale
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
  paste0(sample,"_Organelle0.05Cut.pdf"),
  p,
  width = 7, height = 5,
  units = "in", dpi = 300
)
