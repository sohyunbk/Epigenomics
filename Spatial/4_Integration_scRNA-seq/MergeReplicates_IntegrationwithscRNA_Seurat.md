Merge replicates & integrate with scRNA-seq
================
Sohyun Bang

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: Rcpp

    ## scCustomize v3.0.1
    ## If you find the scCustomize useful please cite.
    ## See 'samuel-marsh.github.io/scCustomize/articles/FAQ.html' for citation info.

    ## qs 0.27.3. Announcement: https://github.com/qsbase/qs/issues/103

Define the Read10X_Image_2 function.. I the error to use “Read10X_Image”
is from the version match error of the spaceranger…

``` r
Read10X_Image_2 <- function(
  image.dir,
  image.name = "tissue_lowres_image.png",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE
) {
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )

  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )

  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(file.path(image.dir, "*tissue_positions*")),
    filter.matrix
  )

  coordinates$imagerow <- as.numeric(coordinates$imagerow)
  coordinates$imagecol <- as.numeric(coordinates$imagecol)

  fov <- CreateFOV(
    coordinates[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = Key(slice, quiet = TRUE)
  )

  visium.fov <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = image,
    scale.factors = scale.factors
  )

  return(visium.fov)
}
```

Load 10X info - image coordinates. Here we do not use Load10X_Spatial…
because of version mismatch issue to read “png” file

``` r
#expr_matrix <- Read10X_h5("./Slide1_WT_Replicate1/filtered_feature_bc_matrix.h5")
features <- read.table("./features.tsv",
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)

CreateSpatialSeuratObject <- function(
  matrix_h5_path,
  image_dir,
  slice = "slice1",
  assay = "Spatial",
  image_function = Read10X_Image_2  # pass your custom image loader
) {
  counts <- Read10X_h5(matrix_h5_path)
  rownames(counts) <- features$V1
  counts <- counts[rownames(counts) != "", ]                    
  counts <- counts[!duplicated(rownames(counts)), ]             
  seurat_obj <- CreateSeuratObject(counts = counts, assay = assay)
  spatial_img <- image_function(image.dir = image_dir, slice = slice)
  DefaultAssay(spatial_img) <- assay
  seurat_obj[[slice]] <- spatial_img
  return(seurat_obj)
  
  }
re1_spatial <- CreateSpatialSeuratObject(
  matrix_h5_path = "./Slide1_WT_Replicate1/filtered_feature_bc_matrix.h5",
  image_dir = "./Slide1_WT_Replicate1/spatial/",
  slice = "Slide1")
```

    ## Warning: Adding image with unordered cells

``` r
re1_spatial$orig.ident <- "re1"

re2_spatial <- CreateSpatialSeuratObject(
  matrix_h5_path = "./Slide2_WT_Replicate2/filtered_feature_bc_matrix.h5",
  image_dir = "./Slide2_WT_Replicate2/spatial/",
  slice = "Slide2")
```

    ## Warning: Adding image with unordered cells

``` r
re2_spatial$orig.ident <- "re2"

SpatialFeaturePlot(re1_spatial, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right") 
```

![](output/fig-Load%20gene%20Info%20data%20&%20QC%20-%20Feature%20and%20count%20number-1.png)<!-- -->

``` r
SpatialFeaturePlot(re2_spatial, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right") 
```

![](output/fig-Load%20gene%20Info%20data%20&%20QC%20-%20Feature%20and%20count%20number-2.png)<!-- -->
Normalization & draw UMAP

``` r
suppressWarnings({
  re1_filt <- SCTransform(re1_spatial, assay = "Spatial", verbose = FALSE)
  re2_filt <- SCTransform(re2_spatial, assay = "Spatial", verbose = FALSE)
})

Merge <- merge(re1_filt, y = re2_filt, add.cell.ids = c("Re1", "Re2"), project = 'maize_visium')
```

    ## Warning: Adding image with unordered cells
    ## Adding image with unordered cells

``` r
features <- SelectIntegrationFeatures(object.list = c(re1_filt, re2_filt), nfeatures = 3000)
VariableFeatures(Merge) <- features

## Dimensionality reduction and clustering
Merge <- Merge %>% 
  RunPCA(assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE, resolution = 0.50) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 12:02:15 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 12:02:15 Read 1583 rows and found 30 numeric columns

    ## 12:02:15 Using Annoy for neighbor search, n_neighbors = 30

    ## 12:02:15 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 12:02:15 Writing NN index file to temp file /var/folders/09/prp0qpds29g8l2p5ts6x3j_m0000gr/T//RtmpIFVW50/file143102227069b
    ## 12:02:15 Searching Annoy index using 1 thread, search_k = 3000
    ## 12:02:15 Annoy recall = 100%
    ## 12:02:15 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 12:02:16 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 12:02:16 Commencing optimization for 500 epochs, with 63720 positive edges
    ## 12:02:16 Using rng type: pcg
    ## 12:02:17 Optimization finished

``` r
DimPlot(Merge, reduction = "umap")
```

![](output/fig-Normalization%20&%20draw%20UMAP-1.png)<!-- -->

``` r
DimPlot(Merge, reduction = "umap", group.by = c("orig.ident"))
```

![](output/fig-Normalization%20&%20draw%20UMAP-2.png)<!-- -->

Find anchors and prediction scores between snRNA-seq and spatial RNA-seq

``` r
snrna_seurat <- readRDS("obj_afterHarmony_WTRe1andRe2_UMI1000.rds")

query_snRNA <- SCTransform(snrna_seurat, verbose = FALSE)

anchors <- FindTransferAnchors(reference = Merge, query = query_snRNA,normalization.method = "SCT")
```

    ## [1] "Given reference assay has multiple sct models, selecting model with most cells for finding transfer anchors"

    ## Selected the SCT model fitted on the most cells.

    ## Performing PCA on the provided reference using 1654 features as input.

    ## Projecting cell embeddings

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 802 anchors

``` r
predictions <- TransferData(
  anchorset = anchors,
  refdata = Merge$seurat_clusters,  # or whatever your annotation column is
  dims = 1:30
)
```

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Predicting cell labels
