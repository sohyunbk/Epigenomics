library(data.table)
library(fgsea)
library(ggplot2)

data(examplePathways)
data(exampleRanks)

str(examplePathways)
head(exampleRanks)
tail(exampleRanks)

fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)