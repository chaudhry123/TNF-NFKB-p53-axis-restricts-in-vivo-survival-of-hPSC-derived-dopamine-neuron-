
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
| bitops 1.0.7| blob 1.2.4 | bluster 1.4.0 | boot 1.3.28.1|



cachem 1.0.8 <br /> circlize 0.4.15 <br /> cli 3.6.1 <br /> clue
0.3.64 <br /> cluster 2.1.4 <br /> codetools 0.2.19 <br /> colorspace
2.1.0 <br /> compiler 4.1.1 <br /> ComplexHeatmap 2.10.0 <br /> cowplot
1.1.1 <br /> crayon 1.5.2 <br /> data.table 1.14.8 <br /> datasets 4.1.1
<br /> DBI 1.1.3 <br /> DelayedArray 0.20.0 <br /> DelayedMatrixStats
1.16.0 <br /> deldir 1.0.6 <br /> digest 0.6.31 <br /> doParallel 1.0.17
<br /> dplyr 1.1.2 <br /> dqrng 0.3.0 <br /> DropletUtils 1.14.2 <br />
edgeR 3.36.0 <br /> ellipsis 0.3.2 <br /> EnhancedVolcano 1.12.0 <br />
evaluate 0.21 <br /> extrafont 0.19 <br /> extrafontdb 1.0 <br /> fansi
1.0.4 <br /> farver 2.1.1 <br /> fastmap 1.1.1 <br /> fitdistrplus
1.1.11 <br /> foreach 1.5.2 <br /> formattable 0.2.1 <br /> future
1.32.0 <br /> future.apply 1.11.0 <br /> generics 0.1.3 <br />
GenomeInfoDb 1.30.1 <br /> GenomeInfoDbData 1.2.7 <br /> GenomicRanges
1.46.1 <br /> GetoptLong 1.0.5 <br /> ggalt 0.4.0 <br /> ggbeeswarm
0.7.2 <br /> ggeasy 0.1.4 <br /> ggplot2 3.4.2 <br /> ggrastr 1.0.2
<br /> ggrepel 0.9.3 <br /> ggridges 0.5.4 <br /> GlobalOptions 0.1.2
<br /> globals 0.16.2 <br /> glue 1.6.2 <br /> goftest 1.2.3 <br />
graphics 4.1.1 <br /> grDevices 4.1.1 <br /> grid 4.1.1 <br /> gridExtra
2.3 <br /> gtable 0.3.3 <br /> HDF5Array 1.22.1 <br /> hms 1.1.3 <br />
htmltools 0.5.5 <br /> htmlwidgets 1.6.2 <br /> httpuv 1.6.10 <br />
httr 1.4.6 <br /> ica 1.0.3 <br /> igraph 1.4.2 <br /> IRanges 2.28.0
<br /> irlba 2.3.5.1 <br /> iterators 1.0.14 <br /> jsonlite 1.8.4
<br /> KEGGREST 1.34.0 <br /> KernSmooth 2.23.21 <br /> knitr 1.43
<br /> labeling 0.4.2 <br /> later 1.3.1 <br /> lattice 0.21.8 <br />
lazyeval 0.2.2 <br /> leiden 0.4.3 <br /> lifecycle 1.0.3 <br /> limma
3.50.3 <br /> listenv 0.9.0 <br /> lme4 1.1.33 <br /> lmtest 0.9.40
<br /> locfit 1.5.9.7 <br /> magick 2.7.4 <br /> magrittr 2.0.3 <br />
maps 3.4.1 <br /> MASS 7.3.60 <br /> MAST 1.20.0 <br /> Matrix 1.5.1
<br /> MatrixGenerics 1.6.0 <br /> matrixStats 0.63.0 <br /> memoise
2.0.1 <br /> metapod 1.2.0 <br /> methods 4.1.1 <br /> mime 0.12 <br />
miniUI 0.1.1.1 <br /> minqa 1.2.5 <br /> monocle3 1.3.1 <br /> munsell
0.5.0 <br /> nlme 3.1.162 <br /> nloptr 2.0.3 <br /> org.Hs.eg.db 3.14.0
<br /> parallel 4.1.1 <br /> parallelly 1.36.0 <br /> patchwork 1.1.2
<br /> pbapply 1.7.0 <br /> pheatmap 1.0.12 <br /> pillar 1.9.0 <br />
pkgconfig 2.0.3 <br /> plotly 4.10.2 <br /> plyr 1.8.8 <br /> png 0.1.8
<br /> polyclip 1.10.4 <br /> prettyunits 1.1.1 <br /> progress 1.2.2
<br /> progressr 0.13.0 <br /> proj4 1.0.12 <br /> promises 1.2.0.1
<br /> purrr 1.0.1 <br /> R.methodsS3 1.8.2 <br /> R.oo 1.25.0 <br />
R.utils 2.12.2 <br /> R6 2.5.1 <br /> RANN 2.6.1 <br /> RColorBrewer
1.1.3 <br /> Rcpp 1.0.10 <br /> RcppAnnoy 0.0.20 <br /> RCurl 1.98.1.12
<br /> reshape2 1.4.4 <br /> reticulate 1.28 <br /> rhdf5 2.38.1 <br />
rhdf5filters 1.6.0 <br /> Rhdf5lib 1.16.0 <br /> rjson 0.2.21 <br />
rlang 1.1.1 <br /> rmarkdown 2.22 <br /> ROCR 1.0.11 <br /> RSQLite
2.3.1 <br /> rstudioapi 0.14 <br /> rsvd 1.0.5 <br /> Rtsne 0.16 <br />
Rttf2pt1 1.3.12 <br /> S4Vectors 0.32.4 <br /> ScaledMatrix 1.2.0 <br />
scales 1.2.1 <br /> scater 1.22.0 <br /> scattermore 0.8 <br /> scran
1.22.1 <br /> sctransform 0.3.5 <br /> scuttle 1.4.0 <br /> Seurat 4.3.0
<br /> SeuratObject 4.1.3 <br /> shape 1.4.6 <br /> shiny 1.7.4 <br />
SingleCellExperiment 1.16.0 <br /> SingleR 1.8.1 <br /> sp 1.6.0 <br />
sparseMatrixStats 1.6.0 <br /> spatstat.data 3.0.1 <br />
spatstat.explore 3.1.0 <br /> spatstat.geom 3.2.1 <br /> spatstat.random
3.1.4 <br /> spatstat.sparse 3.0.1 <br /> spatstat.utils 3.0.3 <br />
splines 4.1.1 <br /> statmod 1.5.0 <br /> stats 4.1.1 <br /> stats4
4.1.1 <br /> stringi 1.7.12 <br /> stringr 1.5.0 <br />
SummarizedExperiment 1.24.0 <br /> survival 3.5.5 <br /> tensor 1.5
<br /> terra 1.7.29 <br /> tibble 3.2.1 <br /> tidyr 1.3.0 <br />
tidyselect 1.2.0 <br /> tools 4.1.1 <br /> utf8 1.2.3 <br /> utils 4.1.1
<br /> uwot 0.1.14 <br /> vctrs 0.6.2 <br /> vipor 0.4.5 <br /> viridis
0.6.3 <br /> viridisLite 0.4.2 <br /> wesanderson 0.3.6 <br /> withr
2.5.0 <br /> xfun 0.39 <br /> xtable 1.8.4 <br /> XVector 0.34.0 <br />
zlibbioc 1.40.0 <br /> zoo 1.8.12 <br />

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
