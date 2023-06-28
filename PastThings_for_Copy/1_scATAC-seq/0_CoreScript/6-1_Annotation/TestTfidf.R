library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
#library(varistran)
library(edgeR)
library(parallel)
library(png)

library("here")
library(devtools)
library(Seurat)
load_all('/home/sb14489/Socrates')
library(dplyr)

meta <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.REF_CELLs.metadata.txt"
geneact <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/GA_A619.txt"
#markers <- "/scratch/sb14489/0.MarkerGene/Ear_All_Combined_Lab_Paper_ReOv.txt"
markers <- "/scratch/sb14489/3.scATAC/0.Data/CellCycle.txt"

pcs <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/Ref/Ref_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt" 

####pcs <- "/scratch/sb14489/3.scATAC_flo/5.Socrates/4_CombineAll_AfterD/Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0_ReduceDSVD_PC100_Var0.REF_CELLs.reduced_dimensions.txt"
#pcs <- "/scratch/sb14489/3.scATAC_flo/5.Socrates/EarStep4_MakeReference_ProjectingMutants/Ref_Tn5Cut10000_BinCut100_MinT0.01_MaxT0.05_ReduceDSVD_PC100_Var50000.REF_CELLs.reduced_dimensions.txt"
threads <- 20
target_cluster <- "LouvainClusters"
plot_each_CT<-"no"


###############
##if we update the markers, we need to re-run this script to obtain the new set of smooth marker file
onlysmooth_marker <- "no" ##yes or no
run_denovo <- "no" ##yes or no
run_known <- "yes" ##yes or no

loadData           <- function(meta,
                               pcs,
                               geneact,
                               markers,
                               target_cluster){
  
  # verbose
  message(" - loading data ...")
  
  # load###
  b <- read.table(meta, header =T)
  #head(b)
  #b <- read.delim(meta)
  #b <- validateCluster(b,target_cluster)
  rownames(b) <- b$cellID
  marker.info <- read.table(markers, header=T)
  rownames(marker.info) <- marker.info$geneID
  h.pcs <- read.table(pcs,header=T)
  head(h.pcs)
  activity <- read.table(geneact,stringsAsFactors = T)
  b <- b[rownames(b) %in% rownames(h.pcs),]
  #head(b)
  h.pcs <- h.pcs[rownames(h.pcs) %in% rownames(b),]
  #head(h.pcs)
  # reformat sparse
  activity <- sparseMatrix(i=as.numeric(activity$V1),
                           j=as.numeric(activity$V2),
                           x=as.numeric(activity$V3),
                           dimnames=list(levels(activity$V1), levels(activity$V2)))
  activity <- activity[,colnames(activity) %in% rownames(b)]
  activity <- activity[,Matrix::colSums(activity>0) > 0]
  activity <- activity[Matrix::rowSums(activity>0) > 0,]
  activity <- activity[,Matrix::colSums(activity>0) > 0]
  
  b <- b[colnames(activity),]
  h.pcs <- h.pcs[colnames(activity),]
  
  # sanity check
  sanity <- 0
  if(sanity > 0){
    print(head(b))
    print(head(marker.info))
    print(head(h.pcs[,1:5]))
    print(head(activity[,1:5]))
  }
  
  # output
  return(list(h.pcs=h.pcs, b=b, activity=activity, marker.info=marker.info))
  
}
str(dat)

###############################################
dat <- loadData(meta, pcs, geneact, markers,target_cluster)
str(dat)
dat$marker.info <- dat$marker.info[c(1:50),]

#runMajorPriori(all.b, all.activity , all.hpcs, marker.info, plot_each_CT,
#               marker_opt_dir,WD,onlysmooth_marker,threads=threads,output=Name)

all.b <- dat$b
all.activity <- dat$activity
rownames(all.activity)[1:100]
all.hpcs <- dat$h.pcs
marker.info <- dat$marker.info

head(all.b)
message("The column name for meta data should be Cluster")
colnames(all.b)[length(colnames(all.b))] 
colnames(all.b)[length(colnames(all.b))] <- "Cluster"
head(all.b)

runMajorPriori     <- function(all.b,
                               all.activity,
                               all.hpcs,
                               marker.info,
                               plot_each_CT,
                               marker_opt_dir,
                               output_dir,
                               onlysmooth_marker,
                               threads=1,
                               output="all",
                               plotHeatmaps=F){
  
  #all.b = metadata
  #all.activity = Cell counts
  # align cell barcodes
  all.activity <- all.activity[,rownames(all.b)]
  all.activity <- all.activity[Matrix::rowSums(all.activity)>0,]
  all.activity <- all.activity[,Matrix::colSums(all.activity)>0]
  head(all.activity)[,c(1:10)]
  #head(all.b)
  ids <- intersect(rownames(all.b), colnames(all.activity))
  ids <- intersect(ids, rownames(all.hpcs))
  #dim(all.b)
  #length(ids)
  all.b <- all.b[ids,]
  all.activity <- all.activity[,ids]
  all.hpcs <- all.hpcs[ids,]
  #head(all.hpcs)
  # normalize per cell activity by cluster average and size factors
  ##updating 071521 check if the normalized activity exits or not
  # load data
  if(file.exists(paste0(output_dir,'/',output,".normalizedActivity_mtx.rds"))){
    activity <- readRDS(paste0(output_dir,'/',output,".normalizedActivity_mtx.rds"))
    #activity <- read.table(paste0(output_dir,"/all.countsActivity.sparse"),stringsAsFactors = T)
  }else{
    results <- normalize.activity(output_dir,all.b,
                                  all.activity,
                                  output=output,
                                  logTransform=F,
                                  scaleP=F,
                                  plotHeatmapRaw=F)
    activity <- results$norm.act
    row.o <- results$row.o
    rm(results)
    gc()
  }
  normalize.activity <- function(output_dir,df,
                                 acts,
                                 n.random=NULL,
                                 output="output",
                                 logTransform=F,
                                 scaleP=F,
                                 plotHeatmapRaw=F)
  df <- all.b
  acts <- all.activity
  head(df)
  head(acts)[,c(1:10)]
  ##########
  ## Here is normalization!!! #####3
  # verbose
  message(" - normalizing libraries ...")
  
  # quick clean
  acts <- acts[Matrix::rowSums(acts>0)>0,]
  acts <- acts[,Matrix::colSums(acts>0)>0]
  #dim(acts)
  # check data frames/matrices
  
  df <- df[rownames(df) %in% colnames(acts),]
  acts <- acts[,colnames(acts) %in% rownames(df)]
  df <- df[colnames(acts),]
  acts.o <- acts
  # if select random
  if(!is.null(n.random)){
    acts <- acts[sample(nrow(acts), n.random),]
  }
  
  # find cluster means
  its <- 0
  message(" - extracting average gene activity across cells per cluster ...")
  
  clust.means <- lapply(unique(df$Cluster), function(x){
    clust.ids <- rownames(subset(df, df$Cluster==x))
    Matrix::rowMeans(acts[,clust.ids])
  })
  clust.ids <- rownames(subset(df, df$Cluster==1))
  tm <- Matrix::rowMeans(acts[,clust.ids])
  length(tm)  
  
  head(df) ## df = meta
  clust.means <- do.call(cbind, clust.means)
  head(clust.means)
  colnames(clust.means) <- unique(df$Cluster)
  
  # plot raw cluster means
  clust.means <- as.matrix(clust.means)
  c.means <- clust.means[,mixedorder(colnames(clust.means))]
  row.o <- apply(c.means, 1, which.max)
  c.means <- c.means[order(row.o, decreasing=F),]
  row.o <- rownames(c.means)
  head(c.means)  
  head(clust.means)
  # if do plot heatmap
  if(plotHeatmapRaw==TRUE){
    pdf(paste0("raw_heatmap.pdf"), width=5, height=6)
    heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
              useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
    dev.off()
  }
  write.table(c.means, file=paste0(output_dir,'/',output,".RawClusterMeans.txt"),
              quote=F, row.names=T, col.names=T, sep="\t")
  
  
  # fit 2-guassian mixture model, return higher threshold
  message(" - fitting mixture model ...")
  thresholds <- apply(clust.means, 2, function(x){
    #x<-clust.means
    #head(x)[c(1:3),]
    #dim(x)
    mod <- Mclust(x, G=2, verbose=F)
    #str(mod)
    #head(mod)[,c(1:10)]
    top <- names(mod$parameters$mean)[which.max(mod$parameters$mean)]
    upper.cells <- x[which(mod$classification==top)]
    val <- quantile(upper.cells[upper.cells>0], c(0.05))
    #val <- mean(upper.cells, na.rm=T)
    if(val == 0){
      val <- min(upper.cells[upper.cells>0])
    }
    return(val)
    #max(mod$parameters$mean)
  })
  str(thresholds)
  head(thresholds)
  # scale cell activity by cluster-specific threshold for expression
  message(" - scaling by cluster averages ...")
  adj.thresh <- thresholds[df$Cluster]
  print(head(adj.thresh, n=10))
  head(acts.o,n=10)
  adj.act <- acts.o %*% Diagonal(x=1/adj.thresh)
  adj.act@x[is.infinite(adj.act@x)] <- 0
  adj.act@x[is.na(adj.act@x)] <- 0
  adj.act@x <- round(adj.act@x, digits=0)
  adj.act <- adj.act[Matrix::rowSums(adj.act>0)>0,]
  adj.act <- adj.act[,Matrix::colSums(adj.act>0)>0]
  
  # un-normalized counts
  ua.out <- as.data.frame(summary(adj.act))
  ua.out$i <- rownames(adj.act)[ua.out$i]
  ua.out$j <- colnames(adj.act)[ua.out$j]
  ua.out <- ua.out[ua.out$x>0,]
  write.table(ua.out, file=paste0(output_dir,'/',output,".countsActivity.sparse"),
              quote=F, row.names=F, col.names=F, sep="\t")
  
  # normalize by size factors
  if(logTransform==F){
    message(" - lib size before size factors = ")
    print(head(Matrix::colSums(adj.act)))
  }
  
  # do normalization
  message(" - estimating normalization factors ...")
  normalizeSF        <- function(adj.act,
                                verbose=F){
    
    # verbose
    if(verbose==T){
      message(" - unnormalized activities:")
      print(head(adj.act[,1:5]))
    }
    
    # get per sample norm factors
    norm.factor <- estimate_sf_sparse(adj.act)
    if(verbose==T){
      message(" - normalization factors:")
      print(head(norm.factor))
    }
    
    # normalized counts
    norm.act <- adj.act %*% Diagonal(x=1/norm.factor)
    colnames(norm.act) <- colnames(adj.act)
    rownames(norm.act) <- rownames(adj.act)
    if(verbose==T){
      message(" - normalized activities:")
      print(head(norm.act[,1:5]))
    }
    
    # output
    return(list(norm.act=norm.act, norm.factor=norm.factor))
  }
  estimate_sf_sparse <- function(counts,
                                 round_exprs=T,
                                 method="mean-geometric-mean-total"){
    if (round_exprs)
      counts <- round(counts)
    
    if(method == 'mean-geometric-mean-total') {
      cell_total <- Matrix::colSums(counts)
      sfs <- cell_total / exp(mean(log(cell_total)))
    }else if(method == 'mean-geometric-mean-log-total') {
      cell_total <- Matrix::colSums(counts)
      sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
    }
    
    sfs[is.na(sfs)] <- 1
    sfs
  }
  results <- normalizeSF(adj.act, verbose=T)
  norm.act <- results$norm.act
  norm.factor <- results$norm.factor
  
  # log transform?
  if(logTransform==T){
    message(" - square-root transforming counts activity ...")
    norm.act <- sqrt(norm.act)
    norm.act@x[is.na(norm.act@x)] <- 0
    norm.act@x[is.infinite(norm.act@x)] <- 0
    message(" - lib size afer square-root transformation = ")
    print(head(Matrix::colSums(norm.act)))
  }
  
  # write size factors to meta data file
  df$size_factors <- norm.factor[rownames(df)]
  head(df) ## df == meta data
  write.table(df, file=paste0(output_dir,'/',output,".size_factors.txt"),
              quote=F, row.names=T, col.names=T, sep="\t")
  
  # verbose
  message(" - lib size after size factors = ")
  print(head(Matrix::colSums(norm.act)))
  
  # remove empty cells/genes
  norm.act <- norm.act[Matrix::rowSums(norm.act>0)>0,]
  norm.act <- norm.act[,Matrix::colSums(norm.act>0)>0]
  message(" - cells = ",ncol(norm.act), " | genes = ", nrow(norm.act))
  print(head(norm.act[,1:10]))
  
  # print output to disk
  ia.out <- as.data.frame(summary(norm.act))
  ia.out$i <- rownames(adj.act)[ia.out$i]
  ia.out$j <- colnames(adj.act)[ia.out$j]
  ia.out <- ia.out[ia.out$x>0,]
  
  
  write.table(ia.out, file=paste0(output_dir,'/',output,".normalizedActivity.sparse"), quote=F, row.names=F,
              col.names=F, sep="\t")
  
  saveRDS(norm.act,paste0(output_dir,'/',output,".normalizedActivity_mtx.rds"))
  
  # return
  message(" - returning normalized matrix ...")
  return(list(norm.act=norm.act, norm.factor=norm.factor, adj.act=adj.act, row.o=row.o))  
  
  results <- list(norm.act=norm.act, norm.factor=norm.factor, adj.act=adj.act, row.o=row.o)
  
  activity <- results$norm.act
  row.o <- results$row.o
  rm(results)
  gc()
  
  
  
  # line-up ids
  
  
  barcode.ids <- intersect(colnames(activity), rownames(all.hpcs))
  barcode.ids <- intersect(barcode.ids, rownames(all.b))
  activity <- activity[,barcode.ids]
  all.hpcs <- all.hpcs[barcode.ids,]
  all.b <- all.b[barcode.ids,]
  
  ##we will smooth all the genes
  # if smooth genes only
  ##since we close this function, we can ignore this and the impute.activity is from all the genes
  if (onlysmooth_marker == 'yes'){
    activity <- activity[rownames(activity) %in% as.character(marker.info$geneID),]
    dim(activity)
  }
  #if(smooth.markers){
  #    activity <- activity[rownames(activity) %in% as.character(marker.info$geneID),]
  #}
  
  # impute gene accessibility scores
  ## Do not run this one
  impute.activity <- smooth.data(activity,
                                 k=15,
                                 step=2,
                                 npcs=ncol(all.hpcs),
                                 df=NULL,
                                 rds=as.data.frame(all.hpcs),
                                 cleanExp=F,
                                 output=output)
  
  #####################################
  ####### smooth data!!
  ######################################
  message(" - imputing gene activity ...")
  smooth.data        <- function(x,
                                 k=25,
                                 step=3,
                                 npcs=19,
                                 cleanExp=F,
                                 df=NULL,
                                 rds=NULL,
                                 output="output")
  # input
  x <- activity
  data.use <- x
  
  # verbose
  if(!is.null(rds)){
    
    if(!is.null(df)){
      message("   * using UMAP manifold for smoothing ...")
      pcs <- df[,c("umap1","umap2")]
    }else{
      message("   * using prior PC space as manifold ...")
      pcs <- rds[colnames(x),c(1:npcs)]
    }
  }else{
    
    # LSI
    message("   * PC manifold set to NULL, running LSI (TFIDF)...")
    x[x>0] <- 1
    tfidf              <- function(bmat,
                                   frequencies=F,
                                   log_scale_tf=T,
                                   scale_factor=100000){
      
      # Use either raw counts or divide by total counts in each cell
      if (frequencies) {
        # "term frequency" method
        tf = t(t(bmat) / Matrix::colSums(bmat))
      } else {
        # "raw count" method
        tf = bmat
      }
      
      # Either TF method can optionally be log scaled
      if (log_scale_tf) {
        if (frequencies) {
          tf@x = log1p(tf@x * scale_factor)
        } else {
          tf@x = log1p(tf@x * 1)
        }
      }
      
      # IDF w/ "inverse document frequency smooth" method
      idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
      
      # TF-IDF
      tf_idf_counts = safe_tfidf(tf, idf)
      rownames(tf_idf_counts) = rownames(bmat)
      colnames(tf_idf_counts) = colnames(bmat)
      return(Matrix(tf_idf_counts, sparse=T))
    }
    safe_tfidf    <- function(tf,
                              idf,
                              block_size=2000e6){
      result = tryCatch({
        result = tf * idf
        result
      }, error = function(e) {
        options(DelayedArray.block.size=block_size)
        DelayedArray:::set_verbose_block_processing(TRUE)
        
        tf = DelayedArray(tf)
        idf = as.matrix(idf)
        
        result = tf * idf
        result
      })
      return(result)
    }
    tf.idf <- tfidf(x)
    head(x)
    # get PCS
    message("   * PC manifold set to NULL, running LSI ...")
    pc <- irlba(t(tf.idf), npcs)
    pcs <- pc$u
    rownames(pcs) <- colnames(x)
    colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
    
    # do l2-norm
    pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
  }
  
  # get KNN
  x1 <- runif(10, 0, 2*pi)
  x2 <- runif(10, 0,3)
  head(x1)
  DATA <- data.frame(x1, x2)
  
  nearest <- nn2(DATA)
  
  message("   * finding knn graph ...")
  head(pcs,n=10)
  dim(pcs)
  #?nn2
  knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
  dim(knn.graph)
  head(knn.graph,n=10)
  j <- as.numeric(x = t(x = knn.graph))
  i <- ((1:length(x = j)) - 1) %/% k + 1
  edgeList = data.frame(i, j, 1)
  A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
  
  # Smooth graph
  ##some questions
  ##1. what's meanning of A
  ##2. A/Matrix::rowSums(A)
  ##3. meanning of A%*%A
  message("   * smoothing graph ...")
  A = A + t(A) ##A now is symmetric
  A = A / Matrix::rowSums(A)
  step.size = step
  if(step.size > 1){
    for(i in 1:step.size){
      message("     ~ step ",i)
      A = A %*% A
    }
  }
  
  # smooth data
  message("   * smoothing activity ...")
  out <- t(A %*% t(data.use))
  colnames(out) <- colnames(x)
  rownames(out) <- rownames(x)
  print(head(out[,1:10]))
  
  # find non-expressed genes
  if(cleanExp==T){
    message("   * finding activity thresholds ...")
    its <- 0
    new <- sparse_apply(out, 2, function(z){
      its <<- its + 1
      mod <- Mclust(z, G=2)
      top.parameter <- which.max(mod$parameters$mean)
      vals <- split(z, mod$classification)[[top.parameter]]
      threshold <- quantile(vals[vals>0], c(0.05))[1]
      if(its %% 100 == 0){message("   * iterated over ",its, " cells (min thresh = ",threshold,") ...")}
      z[z < threshold] <- 0
      return(z)
    }, convert_to_dense=F)
    impute.activity <- Matrix(new, sparse=T)
  }else{
    impute.activity <- out
  }
  
  
  
  
  
  
  
  
  
  
  
  
  ##updating 062221
  if (onlysmooth_marker == 'yes'){
    saveRDS(impute.activity,paste0(output_dir,'/opt_markers_impute.activity.rds'))
  }else{
    saveRDS(impute.activity,paste0(output_dir,'/opt_allgenes_impute.activity.rds'))
  }
  
  # reformat sparse matrix
  #ia <- as.data.frame(summary(impute.activity))
  #ia$i <- rownames(impute.activity)[as.numeric(ia$i)]
  #ia$j <- colnames(impute.activity)[as.numeric(ia$j)]
  #write.table(ia, file=paste0(output,".smoothed.sparse"), quote=F, row.names=F, col.names=F, sep="\t")
  
  # ranges
  marker.impact <- impute.activity[rownames(impute.activity) %in% as.character(marker.info$geneID),]
  marker.ranges <- lapply(rownames(marker.impact), function(x){
    x[is.na(x)] <- 0
    min.x <- min(marker.impact[x,], na.rm=T)
    max.x <- max(marker.impact[x,], na.rm=T)
    rng <- c(min.x,max.x)
    rng[is.na(rng)] <- 0
    rng[1] <- rng[1] - 0.1
    rng[2] <- rng[2] + (0.05*rng[2])
    return(rng)
  })
  names(marker.ranges) <- rownames(marker.impact)
  
  ##updating 062221
  if (plot_each_CT == 'yes'){
    ##plot the markers in each folder
    plot.act.scores_eachCT(marker_opt_dir,all.b, acts=impute.activity,
                           info=marker.info,
                           logT=F,
                           lim=0.999,
                           marker.dist=NULL,
                           outname=paste0("combined.",output,".impute.known.Markers.png"))
    
  }else{
    # plot all
    ##do not add the output_dir to the outname
    plot.act.scores(all.b, output_dir,acts=activity,
                    info=marker.info,
                    logT=T,
                    lim=0.99,
                    marker.dist=NULL,
                    outname=paste0("combined.",output,".normalized.known.Markers.png"))
    
    plot.act.scores(all.b, output_dir,acts=impute.activity,
                    info=marker.info,
                    logT=F,
                    lim=0.999,
                    marker.dist=NULL,
                    outname=paste0("combined.",output,".impute.known.Markers.png"))
    
  }
  # return
  return(list(b=all.b, activity=activity, impute.activity=impute.activity))