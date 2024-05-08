
# TNF NFKB p53 axis restricts in vivo survival of hPSC derived dopamine neuron

This repository contains the code and scripts used for the analysis of
CRISPR screen, bulk RNAseq and scRNAseq generated in:

[Kim TW. et al., TNF NFKB p53 axis restricts in vivo survival of hPSC
derived dopamine neuron, bioRxiv,
2023](https://doi.org/10.1101/2023.03.29.534819) (<font size="2">
*Update link with peer-reviewed publication* </font>)

**Abstract**

Ongoing, first-in-human clinical trials illustrate the feasibility and
translational potential of human pluripotent stem cell (hPSC)-based cell
therapies in Parkinson’s disease (PD). However, a major unresolved
challenge in the field is the extensive cell death following
transplantation with \<10% of grafted dopamine neurons surviving. Here,
we performed a pooled CRISPR/Cas9 screen to enhance survival of
postmitotic dopamine neurons in vivo. We identified p53-mediated
apoptotic cell death as major contributor to dopamine neuron loss and
uncovered a causal link of TNFa-NFκB signaling in limiting cell
survival. As a translationally applicable strategy to purify postmitotic
dopamine neurons, we performed a cell surface marker screen that enabled
purification without the need for genetic reporters. Combining cell
sorting with adalimumab pretreatment, a clinically approved and widely
used TNFa inhibitor, enabled efficient engraftment of postmitotic
dopamine neurons leading to extensive re-innervation and functional
recovery in a preclinical PD mouse model. Thus, transient TNFa
inhibition presents a clinically relevant strategy to enhance survival
and enable engraftment of postmitotic human PSC-derived dopamine neurons
in PD.

**CONTENT**:

-   [Bulk RNA-seq](#Bulk-RNA-seq-data-processing)
-   [Single Cell](#scRNA-seq-data-processing)
-   [Data for download](#data-for-download)

Scripts were written by Fayzan Chaudhry and Hyein Cho.

Samples were processed at the the sequencing facilities at MSKCC and
Weill Cornell Medicine.

Please email us with questions related to the code.

Fayzan Chaudhry
[ffc4001\@med.cornell.edu](mailto:ffc4001@med.cornell.edu){.email}

Hyein Cho [choh\@mskcc.org](mailto:choh@mskcc.org){.email}

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png) \## Bulk RNA-seq data processing

We performed bulk RNAseq analysis to compare gene expression of grafted
NURR1+ neurons, re-isolated at 1 dpt from the mouse brain (d1 graft),
versus that of the matched FACS-purified neurons isolated immediately
prior to transplantation (d0) or analyzed at 1 day of in vitro culture
(d1 culture). Principal component analysis (PCA) and dendrogram analyses
demonstrated that grafted dopamine neurons exhibited the most distinct
transcriptional pattern compared to either sorted or in vitro cultured
dopamine neurons.

Sequences were aligned to *hg19* reference using `STAR` aligner.
Two-pass alignment was used to increase sensitivity. Reads were
quantified and differential expressions were calculated using `DESeq2`.
Functional enrichment analyses were done on the hallmarks of cancer gene
dataset and gene ontology.

## CRISPR

CRISPR screen sequences were quantified by `MAGeCK`, and the differences
and the correlations were generated from the quantifications.

## scRNA-seq data processing {#scrna-seq-data-processing}

Single cell analysis was performed in accordance with the guidelines and
suggestions discussed in the excellent online book *Amezquita R. et al.*
[Amezquita R. et al., Nature Methods, 2019, "Orchestrating single-cell
analysis"](https://osca.bioconductor.org/) describing single-cell
analyses carried out with packages of the Bioconductor environment.

### Processing and analysis

-   The single-cell experiment underwent meticulous processing steps to
    unveil the intricacies of cellular dynamics. Initially, 10X Chromium
    Single Cell 3' v3 technology orchestrated the capture of cellular
    transcripts. These transcripts were then aligned to the human GRCh38
    genome using the venerable `CellRanger 5.0.0`, based on
    `GENCODE v32/Ensembl98` annotations.

-   Following alignment, a stringent filtering process ensued. Cells
    were filtered based on the following criteria:

    -   A minimum of 1000 Unique Molecular Identifier (UMI) counts.
    -   Gene counts ranging from 500 to 7000.
    -   Mitochondrial gene percentage below 25%.

-   Normalization was then achieved through deconvolution using
    `scran V1.22.1`. Cell cycle-related gene expression signature was
    removed using `Seurat V.4.1`.

-   For subsequent analysis 2000 most highly variable genes were
    selected and top 50 principal components were used from the cell
    cycle-regressed matrix.

-   Shared nearest neighbors clustering was constructed using the
    `buildSNNGraph` function in R's `scran`, with k = 40 and cellular
    clusters defined by `walktrap` algorithm implemented in
    `igraph V.1.3.5`.

-   Differential gene expression was elucidated through the MAST
    algorithm in the `Seurat` package.

-   Temporal dynamics were explored through pseudotime analysis via
    `Monocle` in R, while cellular trajectories were delineated using
    `scVelo` in Python. Cluster annotation identified cellular
    identities using `clusterProfiler V.4.2.2`.

-   Visualization of differential expression by volcano plots was
    performed using `EnhancedVolcano V.1.12.0`.

## Data Access

The raw (fastq files) and processed (read counts from HTSeq-count and
CellRanger) can be downloaded from GEO
[GSE216365](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216365).

The easiest way to get started is to use the processed data provided
here.

For the CRISPR screen data, refer to the `internal-data`
[directory](https://github.com/chaudhry123/TNF-NFKB-p53-axis-restricts-in-vivo-survival-of-hPSC-derived-dopamine-neuron-/tree/main/CRISPR-screen-processing/internal-data)
in this repository. Follow through the code highlighted in
[code-to-generate-fig-1.Rmd](https://github.com/chaudhry123/TNF-NFKB-p53-axis-restricts-in-vivo-survival-of-hPSC-derived-dopamine-neuron-/blob/main/CRISPR-screen-processing/code-to-generate-fig-1.Rmd).

For the bulk RNA-seq data, refer to the `internal-data`
[directory](https://github.com/chaudhry123/TNF-NFKB-p53-axis-restricts-in-vivo-survival-of-hPSC-derived-dopamine-neuron-/tree/main/Bulk-RNAseq-processing/internal-data)
in this repository. Follow through the code highlighted in
[bulk-rnaseq-run-through.sh](https://github.com/chaudhry123/TNF-NFKB-p53-axis-restricts-in-vivo-survival-of-hPSC-derived-dopamine-neuron-/blob/main/Bulk-RNAseq-processing/bulk-rnaseq-run-through.sh).

For the single-cell data, some of the data can be downloaded from Box
<https://mskcc.box.com/s/wn5uvwxu2xm4hw219mo0id3r9nkyyprx>

### Software Packages

Example of how to convert a long list of packages to a more readable
table. PLEASE COMPLETE FORMATTING THIS TABLE

|                      |                      |                     |
|----------------------|----------------------|---------------------|
| abind 1.4.5          | AnnotationDbi 1.56.2 | ape 5.7.1           |
| ash 1.0.15           | base 4.1.1           | beachmat 2.10.0     |
| beeswarm 0.4.0       | Biobase 2.54.0       | BiocGenerics 0.40.0 |
| BiocNeighbors 1.12.0 | BiocParallel 1.28.3  | BiocSingular 1.10.0 |
| Biostrings 2.62.0    | bit 4.0.5            | bit64 4.0.5         |
| bitops 1.0.7         | blob 1.2.4           | bluster 1.4.0       |
| boot 1.3.28.1        |cachem 1.0.8          | circlize 0.4.15     | 
| cli 3.6.1            | clue 0.3.64          | cluster 2.1.4       | 
| codetools 0.2.19     | colorspace 2.1.0     | compiler 4.1.1      | 
| ComplexHeatmap 2.10.0| cowplot 1.1.1        | crayon 1.5.2        | 
| data.table 1.14.8    | datasets 4.1.1       | DBI 1.1.3           | 
| DelayedArray 0.20.0  | DelayedMatrixStats 1.16.0| deldir 1.0.6    | 
| digest 0.6.31        | doParallel 1.0.17    | dplyr 1.1.2         | 
| dqrng 0.3.0          | DropletUtils 1.14.2  | edgeR 3.36.0        | 
| ellipsis 0.3.2       | EnhancedVolcano 1.12.0| evaluate 0.21      | 
| extrafont 0.19       | extrafontdb 1.0        | fansi 1.0.4       | 
| farver 2.1.1         | fastmap 1.1.1        | fitdistrplus 1.1.11 |
| foreach 1.5.2        | formattable 0.2.1    | future 1.32.0       | 
| future.apply 1.11.0  | generics 0.1.3       | GenomeInfoDb 1.30.1 |
| GenomeInfoDbData 1.2.7| GenomicRanges 1.46.1 | GetoptLong 1.0.5   |
| ggalt 0.4.0          | ggbeeswarm 0.7.2        | ggeasy 0.1.4        | 
| ggplot2 3.4.2        | ggrastr 1.0.2        | ggrepel 0.9.3        | 
| ggridges 0.5.4       | GlobalOptions 0.1.2  | globals 0.16.2        |
| glue 1.6.2           | goftest 1.2.3  | graphics 4.1.1             |
| grDevices 4.1.1      | grid 4.1.1        | gridExtra 2.3         | 
| gtable 0.3.3        | HDF5Array 1.22.1        | hms 1.1.3  |
| htmltools 0.5.5        | htmlwidgets 1.6.2        | httpuv 1.6.10  |
| httr 1.4.6        | ica 1.0.3        | igraph 1.4.2        |
| IRanges 2.28.0| irlba 2.3.5.1        | iterators 1.0.14      
| jsonlite 1.8.4   | KEGGREST 1.34.0        | KernSmooth 2.23.21       
| knitr 1.43 | labeling 0.4.2        | later 1.3.1        
| lattice 0.21.8  | lazyeval 0.2.2        | leiden 0.4.3        | 
| lifecycle 1.0.3        | limma 3.50.3        | listenv 0.9.0        |
| lme4 1.1.33        | lmtest 0.9.40 | locfit 1.5.9.7        | 
| magick 2.7.4        | magrittr 2.0.3  | maps 3.4.1        | 
| MASS 7.3.60        | MAST 1.20.0        | Matrix 1.5.1| 
| MatrixGenerics 1.6.0        | matrixStats 0.63.0        | memoise2.0.1        |
| metapod 1.2.0        | methods 4.1.1        | mime 0.12  |
| miniUI 0.1.1.1        | minqa 1.2.5        | monocle3 1.3.1        | 
| munsell 0.5.0        | nlme 3.1.162        | nloptr 2.0.3        | 
| org.Hs.eg.db 3.14.0 | parallel 4.1.1        | parallelly 1.36.0        | 
| patchwork 1.1.2 | pbapply 1.7.0        | pheatmap 1.0.12        | 
| pillar 1.9.0  | pkgconfig 2.0.3        | plotly 4.10.2        |
| plyr 1.8.8        | png 0.1.8 | polyclip 1.10.4        |
| prettyunits 1.1.1        | progress 1.2.2 | progressr 0.13.0        |
| proj4 1.0.12        | promises 1.2.0.1 | purrr 1.0.1        | 
| R.methodsS3 1.8.2        | R.oo 1.25.0  | R.utils 2.12.2        |
| R6 2.5.1        | RANN 2.6.1        | RColorBrewer 1.1.3        |
| Rcpp 1.0.10        | RcppAnnoy 0.0.20        | RCurl 1.98.1.12 | 
| reshape2 1.4.4        | reticulate 1.28        | rhdf5 2.38.1  |
| rhdf5filters 1.6.0        | Rhdf5lib 1.16.0        | rjson 0.2.21  |
| rlang 1.1.1        | rmarkdown 2.22        | ROCR 1.0.11        | 
| RSQLite 2.3.1        | rstudioapi 0.14        | rsvd 1.0.5        |
| Rtsne 0.16  | Rttf2pt1 1.3.12        | S4Vectors 0.32.4        |
| ScaledMatrix 1.2.0  | scales 1.2.1        | scater 1.22.0        |
| scattermore 0.8        | scran 1.22.1        | sctransform 0.3.5  |     
| scuttle 1.4.0        | Seurat 4.3.0 | SeuratObject 4.1.3        |
| shape 1.4.6        | shiny 1.7.4  | SingleCellExperiment 1.16.0  |
| SingleR 1.8.1        | sp 1.6.0  | sparseMatrixStats 1.6.0        | 
| spatstat.data 3.0.1  | spatstat.explore 3.1.0        | spatstat.geom 3.2.1        | 
| spatstat.random 3.1.4        | spatstat.sparse 3.0.1        | spatstat.utils 3.0.3  |
| splines 4.1.1        | statmod 1.5.0        | stats 4.1.1        |
| stats4 4.1.1        | stringi 1.7.12        | stringr 1.5.0  |
| SummarizedExperiment 1.24.0        | survival 3.5.5        | tensor 1.5|
       | terra 1.7.29        | tibble 3.2.1        | tidyr 1.3.0  |
| tidyselect 1.2.0        | tools 4.1.1        | utf8 1.2.3        |
| utils 4.1.1 | uwot 0.1.14        | vctrs 0.6.2        
| vipor 0.4.5        | viridis 0.6.3        | viridisLite 0.4.2        |
| wesanderson 0.3.6        | withr 2.5.0        | xfun 0.39        | 
| xtable 1.8.4        | XVector 0.34.0  | zlibbioc 1.40.0        | 
| zoo 1.8.12  | |  | 


## Python

| Package                       | Version            |
|:------------------------------|:-------------------|
| anndata                       | 0.9.2 <br />       |
| anyio                         | 3.6.2 <br />       |
| appnope                       | 0.1.3 <br />       |
| argon2-cffi                   | 21.3.0 <br />      |
| argon2-cffi-bindings          | 21.2.0 <br />      |
| asttokens                     | 2.2.1 <br />       |
| attrs                         | 23.1.0 <br />      |
| backcall                      | 0.2.0 <br />       |
| backports.functools-lru-cache | 1.6.4 <br />       |
| beautifulsoup4                | 4.12.2 <br />      |
| bleach                        | 6.0.0 <br />       |
| certifi                       | 2023.5.7 <br />    |
| cffi                          | 1.15.1 <br />      |
| charset-normalizer            | 3.1.0 <br />       |
| click                         | 8.1.7 <br />       |
| cloudpickle                   | 2.2.1 <br />       |
| comm                          | 0.1.3 <br />       |
| contourpy                     | 1.0.7 <br />       |
| cycler                        | 0.11.0 <br />      |
| debugpy                       | 1.6.7 <br />       |
| decorator                     | 5.1.1 <br />       |
| defusedxml                    | 0.7.1 <br />       |
| entrypoints                   | 0.4 <br />         |
| executing                     | 1.2.0 <br />       |
| fastjsonschema                | 2.16.3 <br />      |
| filelock                      | 3.12.0 <br />      |
| flit_core                     | 3.8.0 <br />       |
| fonttools                     | 4.39.4 <br />      |
| get-annotations               | 0.1.2 <br />       |
| h5py                          | 3.11.0 <br />      |
| idna                          | 3.4 <br />         |
| igraph                        | 0.11.4 <br />      |
| importlib-metadata            | 6.6.0 <br />       |
| importlib-resources           | 5.12.0 <br />      |
| ipykernel                     | 6.23.0 <br />      |
| ipython                       | 8.4.0 <br />       |
| ipython-genutils              | 0.2.0 <br />       |
| ipywidgets                    | 8.0.6 <br />       |
| jedi                          | 0.18.2 <br />      |
| Jinja2                        | 3.1.2 <br />       |
| joblib                        | 1.2.0 <br />       |
| jsonschema                    | 4.17.3 <br />      |
| jupyter                       | 1.0.0 <br />       |
| jupyter_client                | 8.2.0 <br />       |
| jupyter-console               | 6.6.3 <br />       |
| jupyter_core                  | 5.3.0 <br />       |
| jupyter-events                | 0.6.3 <br />       |
| jupyter_server                | 2.5.0 <br />       |
| jupyter_server_terminals      | 0.4.4 <br />       |
| jupyterlab-pygments           | 0.2.2 <br />       |
| jupyterlab-widgets            | 3.0.7 <br />       |
| kiwisolver                    | 1.4.4 <br />       |
| llvmlite                      | 0.36.0 <br />      |
| loompy                        | 3.0.7 <br />       |
| MarkupSafe                    | 2.1.2 <br />       |
| matplotlib                    | 3.7.1 <br />       |
| matplotlib-inline             | 0.1.6 <br />       |
| mistune                       | 2.0.5 <br />       |
| mpmath                        | 1.3.0 <br />       |
| natsort                       | 8.4.0 <br />       |
| nbclassic                     | 1.0.0 <br />       |
| nbclient                      | 0.7.4 <br />       |
| nbconvert                     | 7.4.0 <br />       |
| nbformat                      | 5.8.0 <br />       |
| nest-asyncio                  | 1.5.6 <br />       |
| networkx                      | 3.1 <br />         |
| notebook                      | 6.5.4 <br />       |
| notebook_shim                 | 0.2.3 <br />       |
| numba                         | 0.53.0 <br />      |
| numpy                         | 1.23.0 <br />      |
| numpy-groupies                | 0.9.22 <br />      |
| packaging                     | 23.1 <br />        |
| pandas                        | 2.0.3 <br />       |
| pandocfilters                 | 1.5.0 <br />       |
| parso                         | 0.8.3 <br />       |
| patsy                         | 0.5.6 <br />       |
| pexpect                       | 4.8.0 <br />       |
| pickleshare                   | 0.7.5 <br />       |
| Pillow                        | 9.5.0 <br />       |
| pip                           | 23.1.2 <br />      |
| pkgutil_resolve_name          | 1.3.10 <br />      |
| platformdirs                  | 3.5.0 <br />       |
| ply                           | 3.11 <br />        |
| prometheus-client             | 0.16.0 <br />      |
| prompt-toolkit                | 3.0.38 <br />      |
| psutil                        | 5.9.5 <br />       |
| ptyprocess                    | 0.7.0 <br />       |
| pure-eval                     | 0.2.2 <br />       |
| pycparser                     | 2.21 <br />        |
| Pygments                      | 2.15.1 <br />      |
| pynndescent                   | 0.5.12 <br />      |
| pyobjc-core                   | 9.1.1 <br />       |
| pyobjc-framework-Cocoa        | 9.1.1 <br />       |
| pyparsing                     | 3.0.9 <br />       |
| PyQt5                         | 5.15.7 <br />      |
| PyQt5-sip                     | 12.11.0 <br />     |
| pyreadr                       | 0.4.7 <br />       |
| pyrsistent                    | 0.19.3 <br />      |
| python-dateutil               | 2.8.2 <br />       |
| python-json-logger            | 2.0.7 <br />       |
| pytz                          | 2023.3 <br />      |
| PyYAML                        | 6.0 <br />         |
| pyzmq                         | 25.0.2 <br />      |
| qtconsole                     | 5.4.3 <br />       |
| QtPy                          | 2.3.1 <br />       |
| requests                      | 2.30.0 <br />      |
| rfc3339-validator             | 0.1.4 <br />       |
| rfc3986-validator             | 0.1.1 <br />       |
| scanpy                        | 1.9.8 <br />       |
| scikit-learn                  | 1.2.2 <br />       |
| scipy                         | 1.10.1 <br />      |
| scvelo                        | 0.3.2 <br />       |
| seaborn                       | 0.13.2 <br />      |
| Send2Trash                    | 1.8.2 <br />       |
| session-info                  | 1.0.0 <br />       |
| setuptools                    | 67.7.2 <br />      |
| shap                          | 0.41.0 <br />      |
| sip                           | 6.7.9 <br />       |
| six                           | 1.16.0 <br />      |
| sklearn                       | 0.0.post5 <br />   |
| slicer                        | 0.0.7 <br />       |
| sniffio                       | 1.3.0 <br />       |
| soupsieve                     | 2.3.2.post1 <br /> |
| stack-data                    | 0.6.2 <br />       |
| statsmodels                   | 0.14.1 <br />      |
| stdlib-list                   | 0.10.0 <br />      |
| sympy                         | 1.12 <br />        |
| terminado                     | 0.17.1 <br />      |
| texttable                     | 1.7.0 <br />       |
| threadpoolctl                 | 3.1.0 <br />       |
| tinycss2                      | 1.2.1 <br />       |
| toml                          | 0.10.2 <br />      |
| tomli                         | 2.0.1 <br />       |
| torch                         | 2.0.1 <br />       |
| torchaudio                    | 2.0.2 <br />       |
| torchvision                   | 0.15.2 <br />      |
| tornado                       | 6.3 <br />         |
| tqdm                          | 4.65.0 <br />      |
| traitlets                     | 5.9.0 <br />       |
| typing_extensions             | 4.5.0 <br />       |
| tzdata                        | 2023.3 <br />      |
| umap-learn                    | 0.5.6 <br />       |
| urllib3                       | 2.0.2 <br />       |
| wcwidth                       | 0.2.6 <br />       |
| webencodings                  | 0.5.1 <br />       |
| websocket-client              | 1.5.1 <br />       |
| wheel                         | 0.40.0 <br />      |
| widgetsnbextension            | 4.0.7 <br />       |
| xgboost                       | 1.7.6 <br />       |
| zipp                          | 3.15.0 <br />      |
