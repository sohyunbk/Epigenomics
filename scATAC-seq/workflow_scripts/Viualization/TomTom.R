#https://www.bioconductor.org/packages/release/bioc/vignettes/memes/inst/doc/core_tomtom.html
#conda activate  /home/sb14489/.conda/envs/JASPAR_act
library(memes)
library(magrittr)
options(meme_db = system.file("extdata/flyFactorSurvey_cleaned.meme", 
                              package = "memes", mustWork = TRUE))
