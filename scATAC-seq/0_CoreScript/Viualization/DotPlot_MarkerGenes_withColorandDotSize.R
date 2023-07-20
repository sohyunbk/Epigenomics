# Assuming you have already installed the necessary packages
library(ggplot2)

# Sample data (replace with your actual data)
set.seed(123) # For reproducibility

## In RNA-seq, the dot size is persent expression and Color is average expression.
## In ATAC-seq, the dot size is persent acc and Color is  average acc or zscore.

##################################################
### Make input:  gene - cell_type - Acc percent - zscore.
###################################################

#### 1) Bring zscore file




# Sample gene expression data
genes <- rep(LETTERS[1:10], each = 5)
cell_types <- rep(c("Type1", "Type2", "Type3", "Type4", "Type5"), times = 10)
expression <- rnorm(50, mean = 5, sd = 2)
gene_data <- data.frame(gene = genes, cell_type = cell_types, expression = expression)

# Sample other data
other_cell_types <- c("Type1", "Type2", "Type3", "Type4", "Type5")
other_variable <- rnorm(5, mean = 10, sd = 3)
other_data <- data.frame(cell_type = other_cell_types, other_variable = other_variable)

# Merge both data frames (assuming cell types match in both)
merged_data <- merge(gene_data, other_data, by = "cell_type")

# Plot using ggplot2
# Plot using ggplot2
ggplot(merged_data, aes(x = gene, y = cell_type, 
                        size = expression, color = expression)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +  # Adjust the range for the size of dots
  labs(title = "Dot Plot with Gene Expression and Other Data",
       x = "Gene",
       y = "Cell Type",
       size = "Expression",
       color = "Expression") +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove the grid
        axis.line = element_line(color = "black"),  # Add x-axis and y-axis lines
        panel.border = element_blank()) +  # Remove the panel border
  
  # Change the color scale to range from blue to white to red
  scale_color_gradient2(low = "blue",
                        mid = "white", high = "red", midpoint = 5)

