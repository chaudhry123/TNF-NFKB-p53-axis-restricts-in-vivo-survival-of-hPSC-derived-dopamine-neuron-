
# Bioinformatic methods

>TNF NFKB p53 axis restricts in vivo survival of hPSC derived dopamine neuron

Ongoing, first-in-human clinical trials illustrate the feasibility and translational potential of human pluripotent stem cell (hPSC)-based cell therapies in Parkinson’s disease (PD). However, a major unresolved challenge in the field is the extensive cell death following transplantation with <10% of grafted dopamine neurons surviving. Here, we performed a pooled CRISPR/Cas9 screen to enhance survival of postmitotic dopamine neurons in vivo. We identified p53-mediated apoptotic cell death as major contributor to dopamine neuron loss and uncovered a causal link of TNFa-NFκB signaling in limiting cell survival. As a translationally applicable strategy to purify postmitotic dopamine neurons, we performed a cell surface marker screen that enabled purification without the need for genetic reporters. Combining cell sorting with adalimumab pretreatment, a clinically approved and widely used TNFa inhibitor, enabled efficient engraftment of postmitotic dopamine neurons leading to extensive re-innervation and functional recovery in a preclinical PD mouse model. Thus, transient TNFa inhibition presents a clinically relevant strategy to enhance survival and enable engraftment of postmitotic human PSC-derived dopamine neurons in PD.

Here, we provide an overview of the data and scripts that we used. 

**CONTENT**:

* [Bulk RNA-seq](#Bulk-RNAseq-processingq)
* [CRISPR](#CRISPR-screen-processing)
* [Single Cell](#Single-cell-RNA-seq-Processing)
* [Data for download](#data-for-download)

Scripts were written by Fayzan Chaudhry and Hyein Cho.
Samples were processed at the the sequencing facilities at MSKCC and Weill Cornell Medicine.
Don't hesitate to get in touch with questions related to the code.

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)
## CRISPR 

### Processing 


## Bulk RNA-seq

### Processing 




## scRNA-seq data processing

In principle, we followed the ideas and workflows discussed by Amezquita et al. in the excellent online book ["Orchestrating single-cell analysis"](https://osca.bioconductor.org/) describing single-cell analyses carried out with packages of the Bioconductor environment.

### Processing and analysis
The samples underwent 10X chromium Single Cell 3' v3 processing. The reads were aligned to human GRCh38 (GENCODE v32/Ensembl 98) using Cell Ranger 5.0.0. The resulting filtered count matrix was further filtered for cells with i) minimum 1000 UMI counts, ii) 500 ≤ gene counts ≤ 7000, iii) and mitochondrial gene percentage of less than 25%. Normalization by deconvolution in scran version 1.22.1 was performed and the signal from the gene expression related to the cell cycle was regressed out as directed by Seurat version 4.1. The default 2000 highly variable genes were selected, and the first 50 principal components were extracted from the cell cycle-regressed matrix. Subsequently, the shared nearest neighbors were calculated from the principal components using buildSNNGraph of R software scran using the k parameter of 40. Seven clusters were identified and using the walktrap algorithm, with the function cluster_walktrap of R implementation of the igraph package version 1.3.5. The uniform manifold approximation and projection (UMAP) was performed. Differential gene expression was performed via the Seurat package using MAST. Pseudotime was conducted with Monocle in R, while velocity was conducted with scVelo in python. Cluster annotation was performed via clusterProfiler package version 4.2.2, and differential expression visualization using EnhancedVolcano version 1.12.0.

## Data for download

The raw data (fastq files, read counts from CellRanger) can be downloaded from GEO (INSERT NUMBER HERE).


>The easiest way to get started is to use the processed data provided here.

For the single-cell data, some of the data can be downloaded from Box in the form of RDS (load into R via `in_data <- readRDS()

## Packages
[1] "abind 1.4.5"
[1] "AnnotationDbi 1.56.2"
[1] "ape 5.7.1"
[1] "ash 1.0.15"
[1] "base 4.1.1"
[1] "beachmat 2.10.0"
[1] "beeswarm 0.4.0"
[1] "Biobase 2.54.0"
[1] "BiocGenerics 0.40.0"
[1] "BiocNeighbors 1.12.0"
[1] "BiocParallel 1.28.3"
[1] "BiocSingular 1.10.0"
[1] "Biostrings 2.62.0"
[1] "bit 4.0.5"
[1] "bit64 4.0.5"
[1] "bitops 1.0.7"
[1] "blob 1.2.4"
[1] "bluster 1.4.0"
[1] "boot 1.3.28.1"
[1] "cachem 1.0.8"
[1] "circlize 0.4.15"
[1] "cli 3.6.1"
[1] "clue 0.3.64"
[1] "cluster 2.1.4"
[1] "codetools 0.2.19"
[1] "colorspace 2.1.0"
[1] "compiler 4.1.1"
[1] "ComplexHeatmap 2.10.0"
[1] "cowplot 1.1.1"
[1] "crayon 1.5.2"
[1] "data.table 1.14.8"
[1] "datasets 4.1.1"
[1] "DBI 1.1.3"
[1] "DelayedArray 0.20.0"
[1] "DelayedMatrixStats 1.16.0"
[1] "deldir 1.0.6"
[1] "digest 0.6.31"
[1] "doParallel 1.0.17"
[1] "dplyr 1.1.2"
[1] "dqrng 0.3.0"
[1] "DropletUtils 1.14.2"
[1] "edgeR 3.36.0"
[1] "ellipsis 0.3.2"
[1] "EnhancedVolcano 1.12.0"
[1] "evaluate 0.21"
[1] "extrafont 0.19"
[1] "extrafontdb 1.0"
[1] "fansi 1.0.4"
[1] "farver 2.1.1"
[1] "fastmap 1.1.1"
[1] "fitdistrplus 1.1.11"
[1] "foreach 1.5.2"
[1] "formattable 0.2.1"
[1] "future 1.32.0"
[1] "future.apply 1.11.0"
[1] "generics 0.1.3"
[1] "GenomeInfoDb 1.30.1"
[1] "GenomeInfoDbData 1.2.7"
[1] "GenomicRanges 1.46.1"
[1] "GetoptLong 1.0.5"
[1] "ggalt 0.4.0"
[1] "ggbeeswarm 0.7.2"
[1] "ggeasy 0.1.4"
[1] "ggplot2 3.4.2"
[1] "ggrastr 1.0.2"
[1] "ggrepel 0.9.3"
[1] "ggridges 0.5.4"
[1] "GlobalOptions 0.1.2"
[1] "globals 0.16.2"
[1] "glue 1.6.2"
[1] "goftest 1.2.3"
[1] "graphics 4.1.1"
[1] "grDevices 4.1.1"
[1] "grid 4.1.1"
[1] "gridExtra 2.3"
[1] "gtable 0.3.3"
[1] "HDF5Array 1.22.1"
[1] "hms 1.1.3"
[1] "htmltools 0.5.5"
[1] "htmlwidgets 1.6.2"
[1] "httpuv 1.6.10"
[1] "httr 1.4.6"
[1] "ica 1.0.3"
[1] "igraph 1.4.2"
[1] "IRanges 2.28.0"
[1] "irlba 2.3.5.1"
[1] "iterators 1.0.14"
[1] "jsonlite 1.8.4"
[1] "KEGGREST 1.34.0"
[1] "KernSmooth 2.23.21"
[1] "knitr 1.43"
[1] "labeling 0.4.2"
[1] "later 1.3.1"
[1] "lattice 0.21.8"
[1] "lazyeval 0.2.2"
[1] "leiden 0.4.3"
[1] "lifecycle 1.0.3"
[1] "limma 3.50.3"
[1] "listenv 0.9.0"
[1] "lme4 1.1.33"
[1] "lmtest 0.9.40"
[1] "locfit 1.5.9.7"
[1] "magick 2.7.4"
[1] "magrittr 2.0.3"
[1] "maps 3.4.1"
[1] "MASS 7.3.60"
[1] "MAST 1.20.0"
[1] "Matrix 1.5.1"
[1] "MatrixGenerics 1.6.0"
[1] "matrixStats 0.63.0"
[1] "memoise 2.0.1"
[1] "metapod 1.2.0"
[1] "methods 4.1.1"
[1] "mime 0.12"
[1] "miniUI 0.1.1.1"
[1] "minqa 1.2.5"
[1] "monocle3 1.3.1"
[1] "munsell 0.5.0"
[1] "nlme 3.1.162"
[1] "nloptr 2.0.3"
[1] "org.Hs.eg.db 3.14.0"
[1] "parallel 4.1.1"
[1] "parallelly 1.36.0"
[1] "patchwork 1.1.2"
[1] "pbapply 1.7.0"
[1] "pheatmap 1.0.12"
[1] "pillar 1.9.0"
[1] "pkgconfig 2.0.3"
[1] "plotly 4.10.2"
[1] "plyr 1.8.8"
[1] "png 0.1.8"
[1] "polyclip 1.10.4"
[1] "prettyunits 1.1.1"
[1] "progress 1.2.2"
[1] "progressr 0.13.0"
[1] "proj4 1.0.12"
[1] "promises 1.2.0.1"
[1] "purrr 1.0.1"
[1] "R.methodsS3 1.8.2"
[1] "R.oo 1.25.0"
[1] "R.utils 2.12.2"
[1] "R6 2.5.1"
[1] "RANN 2.6.1"
[1] "RColorBrewer 1.1.3"
[1] "Rcpp 1.0.10"
[1] "RcppAnnoy 0.0.20"
[1] "RCurl 1.98.1.12"
[1] "reshape2 1.4.4"
[1] "reticulate 1.28"
[1] "rhdf5 2.38.1"
[1] "rhdf5filters 1.6.0"
[1] "Rhdf5lib 1.16.0"
[1] "rjson 0.2.21"
[1] "rlang 1.1.1"
[1] "rmarkdown 2.22"
[1] "ROCR 1.0.11"
[1] "RSQLite 2.3.1"
[1] "rstudioapi 0.14"
[1] "rsvd 1.0.5"
[1] "Rtsne 0.16"
[1] "Rttf2pt1 1.3.12"
[1] "S4Vectors 0.32.4"
[1] "ScaledMatrix 1.2.0"
[1] "scales 1.2.1"
[1] "scater 1.22.0"
[1] "scattermore 0.8"
[1] "scran 1.22.1"
[1] "sctransform 0.3.5"
[1] "scuttle 1.4.0"
[1] "Seurat 4.3.0"
[1] "SeuratObject 4.1.3"
[1] "shape 1.4.6"
[1] "shiny 1.7.4"
[1] "SingleCellExperiment 1.16.0"
[1] "SingleR 1.8.1"
[1] "sp 1.6.0"
[1] "sparseMatrixStats 1.6.0"
[1] "spatstat.data 3.0.1"
[1] "spatstat.explore 3.1.0"
[1] "spatstat.geom 3.2.1"
[1] "spatstat.random 3.1.4"
[1] "spatstat.sparse 3.0.1"
[1] "spatstat.utils 3.0.3"
[1] "splines 4.1.1"
[1] "statmod 1.5.0"
[1] "stats 4.1.1"
[1] "stats4 4.1.1"
[1] "stringi 1.7.12"
[1] "stringr 1.5.0"
[1] "SummarizedExperiment 1.24.0"
[1] "survival 3.5.5"
[1] "tensor 1.5"
[1] "terra 1.7.29"
[1] "tibble 3.2.1"
[1] "tidyr 1.3.0"
[1] "tidyselect 1.2.0"
[1] "tools 4.1.1"
[1] "utf8 1.2.3"
[1] "utils 4.1.1"
[1] "uwot 0.1.14"
[1] "vctrs 0.6.2"
[1] "vipor 0.4.5"
[1] "viridis 0.6.3"
[1] "viridisLite 0.4.2"
[1] "wesanderson 0.3.6"
[1] "withr 2.5.0"
[1] "xfun 0.39"
[1] "xtable 1.8.4"
[1] "XVector 0.34.0"
[1] "zlibbioc 1.40.0"
[1] "zoo 1.8.12"


