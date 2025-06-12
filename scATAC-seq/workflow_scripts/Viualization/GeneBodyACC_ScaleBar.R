library(RColorBrewer)
cols <- colorRampPalette(c("grey80", "grey76", "grey72", brewer.pal(9, "RdPu")[3:9]), bias = 0.5)(100)

pdf("/Users/sohyun/Documents/2.SingleCellATAC/3.Data/6.Annotation/1.MarkerGene/ScaleBar.pdf", width = 8, height = 4)
barplot(rep(1, length(cols)), space = 0, col = cols, border = NA, xlab = "Index", ylab = "Color Value", main = "Color Palette Display")
dev.off()
