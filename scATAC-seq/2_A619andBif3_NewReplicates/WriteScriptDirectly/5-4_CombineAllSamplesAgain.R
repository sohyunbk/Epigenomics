library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library("optparse")
library(rlang)
library(ggplot2)
library("RColorBrewer")

Ex <- function(){
  Dir <- as.character("CombineAll") #NewDir name
  Prefix <- "Tn5Cut1000_Binsize500_MinT0.005_MaxT0.05_PC100"
  WD <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AfterMtMapping/"
}

args <- commandArgs(T)
Dir <- as.character(args[1])
minimumtn5counts <- as.character(args[2])
Minimum_PeakNumberPerCell <- as.character(args[3])
MinT<- as.character(args[4])
MaxT <- as.character(args[5])
SVDorNMF <- as.character(args[6])
NumberOfPC <- as.character(args[7])

#Ex <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619_2/1_A619_2_Tn5Cut1000_Binsize500_Mt0.05.rds")
#Ex2 <- readRDS("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/1_A619_2/1_A619_2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.rds")

# functions --------------------------------------------------------------

# remove high conf cells from raw data set
LoadPreviousRds <- function(Name,Prefix,libName,NameSample){
  File <- paste0(WD,"/",Name,"/",Name,"_",Prefix)
  Pastobj <- readRDS(paste(File,".rds",sep=""))
  Pastobj$meta$SampleName <- NameSample
  Pastobj$meta$library <- libName
  newobj <- list()
  newobj$meta <- Pastobj$meta 
  newobj$counts <- Pastobj$counts 
  return(newobj)
  }

NotIntersectNumber <- function(Avector,Bvector){
  print(length(union(Avector,Bvector))-length(intersect(Avector,Bvector)))}

#############===================================================================
## Combine All together in the begening!!

obj_A619_Re1 <- LoadPreviousRds("A619_Re3",Prefix,"A619_Re3","A619")
obj_A619_Re2 <- LoadPreviousRds("A619_Re4",Prefix,"A619_Re4","A619")
## tfidf was performed by eash sample. I will combine the replicates and then tfidf again
## Extract only meta/counts for merging! ##
#obj_A619_Re1$meta$lib_ID <- "1"
str(obj_A619_Re1)
str(obj_A619_Re2)
##Check if there are some rows in zero sum
#rownames(obj_A619_Re2$counts)[1:4]
#Temp <- obj_A619_Re2$counts[Matrix::rowSums(obj_A619_Re2$counts) > 0,]


obj_bif3_Re1 <- LoadPreviousRds("bif3_Re3",Prefix,"bif3_Re3","bif3")
obj_bif3_Re2 <- LoadPreviousRds("bif3_Re4",Prefix,"bif3_Re4","bif3")
str(obj_bif3_Re1)
str(obj_bif3_Re2)

#rownames(obj_relk1_Re1$counts)


NotIntersectNumber(rownames(obj_A619_Re1$counts),rownames(obj_A619_Re2$counts))
NotIntersectNumber(rownames(obj_bif3_Re1$counts),rownames(obj_bif3_Re2$counts))

SharedFeatures <- Reduce(intersect, list(rownames(obj_A619_Re1$counts),rownames(obj_A619_Re2$counts),
                                      rownames(obj_bif3_Re1$counts),rownames(obj_bif3_Re2$counts)))

length(SharedFeatures)
print(SharedFeatures[1:5])
#print(SharedFeatures)

files <- list(obj_A619_Re1,obj_A619_Re2, 
              obj_bif3_Re1, obj_bif3_Re2)
names(files) <- c("A619_Re3", "A619_Re4",
                  "bif3_Re3","bif3_Re4")

merged.obj <- mergeSocratesRDS(obj.list=files)

#######################
setwd(WD)
## Make directory
if (file.exists(file.path(Dir))){
  setwd(file.path(Dir))
} else {
  dir.create(file.path(Dir))
  setwd(file.path(Dir))
}
getwd()
str(merged.obj)
#####################
str(merged.obj)
head(merged.obj$counts[,c(1:10)])
merged.obj$counts <- merged.obj$counts[rownames(merged.obj$counts) %in% SharedFeatures,]
str(merged.obj)

Name <- paste0("Combined_",Prefix)
saveRDS(merged.obj, file=paste0(Name,".rds"))

