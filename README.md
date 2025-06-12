# Overview of single cell genomics repo
The scrtips are for the [manuscript](https://www.biorxiv.org/content/10.1101/2024.05.13.593957v1) entitled "WUSCHEL-dependent chromatin regulation in maize inflorescence development at single-cell resolution".
They can be used for more than two genotype comparisions in scRNA-seq, spatial RNA-seq, and scATAC-seq data.
Omics analyses such as ChIP-seq and DAP-seq are also integrated.

# Repository Structure 
## Spatial Transcriptome Analysis
`Spatial` dir contains scripts used for spatial transcriptomics analysis (10X Visium) in maize inflorescence tissue. 
We used spaceranger to produce spot level gene expression data with spatial coordinates.

The pipeline is organized by analysis stages:
* `1_QC_MarkerGene/`: Quality control, normalization, and marker gene identification pipeline using `Scanpy`.
  In the begining, QC & MarkerGene finding step was 
* `2_TargetGene_Replicates/`: Scripts for extracting metadata, summarizing gene counts, and performing differential expression tests on selected target genes across replicates. Includes both Python and R scripts.
* `3_SelectCells/`:
domain selection and position-based extraction of gene counts. Includes DE analysis for selected gene sets (e.g., WOX genes), with interactive notebooks and R scripts for plotting and statistics.

## scRNA-seq Analysis
`scRNA-seq` dir has scripts used for scRNA-seq analysis. 
It includes reference preparation, read alignment, quality control, data integration, and cell type annotation. The pipeline supports both `Cell Ranger` and `STARsolo` aligners and uses `Seurat` and `Harmony` for downstream analyses. 
The pipeline is structured with both:
  - Shell (.sh) scripts: For SLURM submission.
  - R (.R) scripts: For data analysis and visualization within R, including quality control, normalization, integration, and annotation.

## scATAC-seq Analysis
`scATAC-seq` dir contains scripts used for scATAC-seq analysis. 
### submission_scripts
* `submission_scripts` is main directory includes SLURM submission scripts for running the pipeline steps on a compute cluster. Each script wraps a specific component in handling resource allocation, environment setup, and job execution.
They are organized numerically to reflect the order of execution.
* `submission_scripts_additionalsamples` directory has several trials and fails. First, we produced more libaraies before so tried to combine all libraries together. Second, to compare the mutant and wild type there are several method we can use and one of them is to use wild type as reference and project mutant cells to the reference. `3_A619andBif3ToBif3Ref` has the reference based mutant projection method.

### workflow_scripts
`workflow_scripts` directory contains core R, Python, and shell scripts that perform specific analyses in a modular and reusable fashion. Each subdirectory represents a distinct step or module of the analysis pipeline.

* `Alignment/`: Scripts for aligning raw reads using Cell Ranger ATAC, and summarizing mapping results.
* `Annotation_Cluster/`: Scripts for annotating clusters, generating UMAP visualizations, and identifying marker genes.
* `Clustering/`: R scripts to perform clustering, UMAP projection, and comparisons across replicates.
* `Correlation/`: Scripts to calculate correlation matrices across replicates, gene bodies, and peak-level accessibility.
* `dACR/`: Main scripts for identifying and analyzing differentially accessible chromatin regions (dACRs).
* `dACRAdditionalAnalysis/`: Downstream analysis of dACRs, including gene association, heatmaps, and statistical summarization.
* `FindCommonACRPos_with500pbLength/`:  Scripts for identifying reproducible ACR positions and processing them to standardized lengths.
* `Mapping_RefiningBam/`: Shell scripts to refine BAM files via deduplication and removal of multimapped reads.
Meme/: Scripts to discover and analyze motifs in ACRs using MEME suite tools.
* `MotifDeviation/`: R scripts to compute motif accessibility deviations using ChromVAR.
* `PeakCalling_byCellTypes/`: Scripts for cell-type-specific peak calling, merging, and classification of ACRs.
* `SocratesStart_QC/`: Scripts to initialize Socrates objects, filter cells, and assess quality control metrics.
* `Viualization/`: R scripts to generate summary figures including dot plots, density plots, pie charts, and heatmaps.

### functions
`functions` directory contains standalone R functions that are sourced by other scripts. 

### Revision_GenomeBiology
`Revision_GenomeBiology` directory contains the additional analysis & visualization scripts. They were added during revision in genome biology journal.
The scripts are written in Rmd to make it easy to check figures and scripts together.

