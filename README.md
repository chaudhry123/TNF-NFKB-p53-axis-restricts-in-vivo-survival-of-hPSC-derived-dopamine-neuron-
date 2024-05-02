
# Bioinformatic methods

>TNF NFKB p53 axis restricts in vivo survival of hPSC derived dopamine neuron

Ongoing, first-in-human clinical trials illustrate the feasibility and translational potential of human pluripotent stem cell (hPSC)-based cell therapies in Parkinson’s disease (PD). However, a major unresolved challenge in the field is the extensive cell death following transplantation with <10% of grafted dopamine neurons surviving. Here, we performed a pooled CRISPR/Cas9 screen to enhance survival of postmitotic dopamine neurons in vivo. We identified p53-mediated apoptotic cell death as major contributor to dopamine neuron loss and uncovered a causal link of TNFa-NFκB signaling in limiting cell survival. As a translationally applicable strategy to purify postmitotic dopamine neurons, we performed a cell surface marker screen that enabled purification without the need for genetic reporters. Combining cell sorting with adalimumab pretreatment, a clinically approved and widely used TNFa inhibitor, enabled efficient engraftment of postmitotic dopamine neurons leading to extensive re-innervation and functional recovery in a preclinical PD mouse model. Thus, transient TNFa inhibition presents a clinically relevant strategy to enhance survival and enable engraftment of postmitotic human PSC-derived dopamine neurons in PD.

Here, we provide an overview of the data and scripts that we used. 

**CONTENT**:

* [Bulk RNA-seq](#Bulk-RNAseq-processingq)
* [Single Cell](#Single-cell-RNA-seq-Processing)
* [Data for download](#data-for-download)

Scripts were written by Fayzan Chaudhry and Hyein Cho.
Samples were processed at the the sequencing facilities at MSKCC and Weill Cornell Medicine.
Don't hesitate to get in touch with questions related to the code.

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)
## Bulk RNA-seq
We performed bulk RNAseq analysis to compare gene expression of grafted NURR1+ neurons, re-isolated at 1 dpt from the mouse brain (d1 graft), versus that of the matched FACS-purified neurons isolated immediately prior to transplantation (d0) or analyzed at 1 day of in vitro culture (d1 culture). Principal component analysis (PCA) and dendrogram analyses demonstrated that grafted dopamine neurons exhibited the most distinct transcriptional pattern compared to either sorted or in vitro cultured dopamine neurons

## CRISPR

## scRNA-seq data processing

In principle, we followed the ideas and workflows discussed by Amezquita et al. in the excellent online book ["Orchestrating single-cell analysis"](https://osca.bioconductor.org/) describing single-cell analyses carried out with packages of the Bioconductor environment.

### Processing and analysis
The samples underwent 10X chromium Single Cell 3' v3 processing. The reads were aligned to human GRCh38 (GENCODE v32/Ensembl 98) using Cell Ranger 5.0.0. The resulting filtered count matrix was further filtered for cells with i) minimum 1000 UMI counts, ii) 500 ≤ gene counts ≤ 7000, iii) and mitochondrial gene percentage of less than 25%. Normalization by deconvolution in scran version 1.22.1 was performed and the signal from the gene expression related to the cell cycle was regressed out as directed by Seurat version 4.1. The default 2000 highly variable genes were selected, and the first 50 principal components were extracted from the cell cycle-regressed matrix. Subsequently, the shared nearest neighbors were calculated from the principal components using buildSNNGraph of R software scran using the k parameter of 40. Seven clusters were identified and using the walktrap algorithm, with the function cluster_walktrap of R implementation of the igraph package version 1.3.5. The uniform manifold approximation and projection (UMAP) was performed. Differential gene expression was performed via the Seurat package using MAST. Pseudotime was conducted with Monocle in R, while velocity was conducted with scVelo in python. Cluster annotation was performed via clusterProfiler package version 4.2.2, and differential expression visualization using EnhancedVolcano version 1.12.0.

## Data for download

The raw data (fastq files, read counts from CellRanger) can be downloaded from GEO GSE217131.


>The easiest way to get started is to use the processed data provided here.

For the single-cell data, some of the data can be downloaded from Box https://mskcc.box.com/s/dm20vn2ww1c1xto87xmzqb9kxilgxpcc

## Packages
[1] "abind 1.4.5" <br />
[1] "AnnotationDbi 1.56.2" <br />
[1] "ape 5.7.1" <br />
[1] "ash 1.0.15" <br />
[1] "base 4.1.1" <br />
[1] "beachmat 2.10.0" <br />
[1] "beeswarm 0.4.0" <br />
[1] "Biobase 2.54.0" <br />
[1] "BiocGenerics 0.40.0" <br />
[1] "BiocNeighbors 1.12.0" <br />
[1] "BiocParallel 1.28.3" <br />
[1] "BiocSingular 1.10.0" <br />
[1] "Biostrings 2.62.0" <br />
[1] "bit 4.0.5" <br />
[1] "bit64 4.0.5" <br />
[1] "bitops 1.0.7" <br />
[1] "blob 1.2.4" <br />
[1] "bluster 1.4.0" <br />
[1] "boot 1.3.28.1" <br />
[1] "cachem 1.0.8" <br />
[1] "circlize 0.4.15" <br />
[1] "cli 3.6.1" <br />
[1] "clue 0.3.64" <br />
[1] "cluster 2.1.4" <br />
[1] "codetools 0.2.19" <br />
[1] "colorspace 2.1.0" <br />
[1] "compiler 4.1.1" <br />
[1] "ComplexHeatmap 2.10.0" <br />
[1] "cowplot 1.1.1" <br />
[1] "crayon 1.5.2" <br />
[1] "data.table 1.14.8" <br />
[1] "datasets 4.1.1" <br />
[1] "DBI 1.1.3" <br />
[1] "DelayedArray 0.20.0" <br />
[1] "DelayedMatrixStats 1.16.0" <br />
[1] "deldir 1.0.6" <br />
[1] "digest 0.6.31" <br />
[1] "doParallel 1.0.17" <br />
[1] "dplyr 1.1.2" <br />
[1] "dqrng 0.3.0" <br />
[1] "DropletUtils 1.14.2" <br />
[1] "edgeR 3.36.0" <br />
[1] "ellipsis 0.3.2" <br />
[1] "EnhancedVolcano 1.12.0" <br />
[1] "evaluate 0.21" <br />
[1] "extrafont 0.19" <br />
[1] "extrafontdb 1.0" <br />
[1] "fansi 1.0.4" <br />
[1] "farver 2.1.1" <br />
[1] "fastmap 1.1.1" <br />
[1] "fitdistrplus 1.1.11" <br />
[1] "foreach 1.5.2" <br />
[1] "formattable 0.2.1" <br />
[1] "future 1.32.0" <br />
[1] "future.apply 1.11.0" <br />
[1] "generics 0.1.3" <br />
[1] "GenomeInfoDb 1.30.1" <br />
[1] "GenomeInfoDbData 1.2.7" <br />
[1] "GenomicRanges 1.46.1" <br />
[1] "GetoptLong 1.0.5" <br />
[1] "ggalt 0.4.0" <br />
[1] "ggbeeswarm 0.7.2" <br />
[1] "ggeasy 0.1.4" <br />
[1] "ggplot2 3.4.2" <br />
[1] "ggrastr 1.0.2" <br />
[1] "ggrepel 0.9.3" <br />
[1] "ggridges 0.5.4" <br />
[1] "GlobalOptions 0.1.2" <br />
[1] "globals 0.16.2" <br />
[1] "glue 1.6.2" <br />
[1] "goftest 1.2.3" <br />
[1] "graphics 4.1.1" <br />
[1] "grDevices 4.1.1" <br />
[1] "grid 4.1.1" <br />
[1] "gridExtra 2.3" <br />
[1] "gtable 0.3.3" <br />
[1] "HDF5Array 1.22.1" <br />
[1] "hms 1.1.3" <br />
[1] "htmltools 0.5.5" <br />
[1] "htmlwidgets 1.6.2" <br />
[1] "httpuv 1.6.10" <br />
[1] "httr 1.4.6" <br />
[1] "ica 1.0.3" <br />
[1] "igraph 1.4.2" <br />
[1] "IRanges 2.28.0" <br />
[1] "irlba 2.3.5.1" <br />
[1] "iterators 1.0.14" <br />
[1] "jsonlite 1.8.4" <br />
[1] "KEGGREST 1.34.0" <br />
[1] "KernSmooth 2.23.21" <br />
[1] "knitr 1.43" <br />
[1] "labeling 0.4.2" <br />
[1] "later 1.3.1" <br />
[1] "lattice 0.21.8" <br />
[1] "lazyeval 0.2.2" <br />
[1] "leiden 0.4.3" <br />
[1] "lifecycle 1.0.3" <br />
[1] "limma 3.50.3" <br />
[1] "listenv 0.9.0" <br />
[1] "lme4 1.1.33" <br />
[1] "lmtest 0.9.40" <br />
[1] "locfit 1.5.9.7" <br />
[1] "magick 2.7.4" <br />
[1] "magrittr 2.0.3" <br />
[1] "maps 3.4.1" <br />
[1] "MASS 7.3.60" <br />
[1] "MAST 1.20.0" <br />
[1] "Matrix 1.5.1" <br />
[1] "MatrixGenerics 1.6.0" <br />
[1] "matrixStats 0.63.0" <br />
[1] "memoise 2.0.1" <br />
[1] "metapod 1.2.0" <br />
[1] "methods 4.1.1" <br />
[1] "mime 0.12" <br />
[1] "miniUI 0.1.1.1" <br />
[1] "minqa 1.2.5" <br />
[1] "monocle3 1.3.1" <br />
[1] "munsell 0.5.0" <br />
[1] "nlme 3.1.162" <br />
[1] "nloptr 2.0.3" <br />
[1] "org.Hs.eg.db 3.14.0" <br />
[1] "parallel 4.1.1" <br />
[1] "parallelly 1.36.0" <br />
[1] "patchwork 1.1.2" <br />
[1] "pbapply 1.7.0" <br />
[1] "pheatmap 1.0.12" <br />
[1] "pillar 1.9.0" <br />
[1] "pkgconfig 2.0.3" <br />
[1] "plotly 4.10.2" <br />
[1] "plyr 1.8.8" <br />
[1] "png 0.1.8" <br />
[1] "polyclip 1.10.4" <br />
[1] "prettyunits 1.1.1" <br />
[1] "progress 1.2.2" <br />
[1] "progressr 0.13.0" <br />
[1] "proj4 1.0.12" <br />
[1] "promises 1.2.0.1" <br />
[1] "purrr 1.0.1" <br />
[1] "R.methodsS3 1.8.2" <br />
[1] "R.oo 1.25.0" <br />
[1] "R.utils 2.12.2" <br />
[1] "R6 2.5.1" <br />
[1] "RANN 2.6.1" <br />
[1] "RColorBrewer 1.1.3" <br />
[1] "Rcpp 1.0.10" <br />
[1] "RcppAnnoy 0.0.20" <br />
[1] "RCurl 1.98.1.12" <br />
[1] "reshape2 1.4.4" <br />
[1] "reticulate 1.28" <br />
[1] "rhdf5 2.38.1" <br />
[1] "rhdf5filters 1.6.0" <br />
[1] "Rhdf5lib 1.16.0" <br />
[1] "rjson 0.2.21" <br />
[1] "rlang 1.1.1" <br />
[1] "rmarkdown 2.22" <br />
[1] "ROCR 1.0.11" <br />
[1] "RSQLite 2.3.1" <br />
[1] "rstudioapi 0.14" <br />
[1] "rsvd 1.0.5" <br />
[1] "Rtsne 0.16" <br />
[1] "Rttf2pt1 1.3.12" <br />
[1] "S4Vectors 0.32.4" <br />
[1] "ScaledMatrix 1.2.0" <br />
[1] "scales 1.2.1" <br />
[1] "scater 1.22.0" <br />
[1] "scattermore 0.8" <br />
[1] "scran 1.22.1" <br />
[1] "sctransform 0.3.5" <br />
[1] "scuttle 1.4.0" <br />
[1] "Seurat 4.3.0" <br />
[1] "SeuratObject 4.1.3" <br />
[1] "shape 1.4.6" <br />
[1] "shiny 1.7.4" <br />
[1] "SingleCellExperiment 1.16.0" <br />
[1] "SingleR 1.8.1" <br />
[1] "sp 1.6.0" <br />
[1] "sparseMatrixStats 1.6.0" <br />
[1] "spatstat.data 3.0.1" <br />
[1] "spatstat.explore 3.1.0" <br />
[1] "spatstat.geom 3.2.1" <br />
[1] "spatstat.random 3.1.4" <br />
[1] "spatstat.sparse 3.0.1" <br />
[1] "spatstat.utils 3.0.3" <br />
[1] "splines 4.1.1" <br />
[1] "statmod 1.5.0" <br />
[1] "stats 4.1.1" <br />
[1] "stats4 4.1.1" <br />
[1] "stringi 1.7.12" <br />
[1] "stringr 1.5.0" <br />
[1] "SummarizedExperiment 1.24.0" <br />
[1] "survival 3.5.5" <br />
[1] "tensor 1.5" <br />
[1] "terra 1.7.29" <br />
[1] "tibble 3.2.1" <br />
[1] "tidyr 1.3.0" <br />
[1] "tidyselect 1.2.0" <br />
[1] "tools 4.1.1" <br />
[1] "utf8 1.2.3" <br />
[1] "utils 4.1.1" <br />
[1] "uwot 0.1.14" <br />
[1] "vctrs 0.6.2" <br />
[1] "vipor 0.4.5" <br />
[1] "viridis 0.6.3" <br />
[1] "viridisLite 0.4.2" <br />
[1] "wesanderson 0.3.6" <br />
[1] "withr 2.5.0" <br />
[1] "xfun 0.39" <br />
[1] "xtable 1.8.4" <br />
[1] "XVector 0.34.0" <br />
[1] "zlibbioc 1.40.0" <br />
[1] "zoo 1.8.12" <br />

## Python
Package                       Version
----------------------------- -----------
anndata                       0.9.2 <br />
anyio                         3.6.2 <br />
appnope                       0.1.3 <br />
argon2-cffi                   21.3.0 <br />
argon2-cffi-bindings          21.2.0 <br />
asttokens                     2.2.1 <br />
attrs                         23.1.0 <br />
backcall                      0.2.0 <br />
backports.functools-lru-cache 1.6.4 <br />
beautifulsoup4                4.12.2 <br />
bleach                        6.0.0 <br />
certifi                       2023.5.7 <br />
cffi                          1.15.1 <br />
charset-normalizer            3.1.0 <br />
click                         8.1.7 <br />
cloudpickle                   2.2.1 <br />
comm                          0.1.3 <br />
contourpy                     1.0.7 <br />
cycler                        0.11.0 <br />
debugpy                       1.6.7 <br />
decorator                     5.1.1 <br />
defusedxml                    0.7.1 <br />
entrypoints                   0.4 <br />
executing                     1.2.0 <br />
fastjsonschema                2.16.3 <br />
filelock                      3.12.0 <br />
flit_core                     3.8.0 <br />
fonttools                     4.39.4 <br />
get-annotations               0.1.2 <br />
h5py                          3.11.0 <br />
idna                          3.4 <br />
igraph                        0.11.4 <br />
importlib-metadata            6.6.0 <br />
importlib-resources           5.12.0 <br />
ipykernel                     6.23.0 <br />
ipython                       8.4.0 <br />
ipython-genutils              0.2.0 <br />
ipywidgets                    8.0.6 <br />
jedi                          0.18.2 <br />
Jinja2                        3.1.2 <br />
joblib                        1.2.0 <br />
jsonschema                    4.17.3 <br />
jupyter                       1.0.0 <br />
jupyter_client                8.2.0 <br />
jupyter-console               6.6.3 <br />
jupyter_core                  5.3.0 <br />
jupyter-events                0.6.3 <br />
jupyter_server                2.5.0 <br />
jupyter_server_terminals      0.4.4 <br />
jupyterlab-pygments           0.2.2 <br />
jupyterlab-widgets            3.0.7 <br />
kiwisolver                    1.4.4 <br />
llvmlite                      0.36.0 <br />
loompy                        3.0.7 <br />
MarkupSafe                    2.1.2 <br />
matplotlib                    3.7.1 <br />
matplotlib-inline             0.1.6 <br />
mistune                       2.0.5 <br />
mpmath                        1.3.0 <br />
natsort                       8.4.0 <br />
nbclassic                     1.0.0 <br />
nbclient                      0.7.4 <br />
nbconvert                     7.4.0 <br />
nbformat                      5.8.0 <br />
nest-asyncio                  1.5.6 <br />
networkx                      3.1 <br />
notebook                      6.5.4 <br />
notebook_shim                 0.2.3 <br />
numba                         0.53.0 <br />
numpy                         1.23.0 <br />
numpy-groupies                0.9.22 <br />
packaging                     23.1 <br />
pandas                        2.0.3 <br />
pandocfilters                 1.5.0 <br />
parso                         0.8.3 <br />
patsy                         0.5.6 <br />
pexpect                       4.8.0 <br />
pickleshare                   0.7.5 <br />
Pillow                        9.5.0 <br />
pip                           23.1.2 <br />
pkgutil_resolve_name          1.3.10 <br />
platformdirs                  3.5.0 <br />
ply                           3.11 <br />
prometheus-client             0.16.0 <br />
prompt-toolkit                3.0.38 <br />
psutil                        5.9.5 <br />
ptyprocess                    0.7.0 <br />
pure-eval                     0.2.2 <br />
pycparser                     2.21 <br />
Pygments                      2.15.1 <br />
pynndescent                   0.5.12 <br />
pyobjc-core                   9.1.1 <br />
pyobjc-framework-Cocoa        9.1.1 <br />
pyparsing                     3.0.9 <br />
PyQt5                         5.15.7 <br />
PyQt5-sip                     12.11.0 <br />
pyreadr                       0.4.7 <br />
pyrsistent                    0.19.3 <br />
python-dateutil               2.8.2 <br />
python-json-logger            2.0.7 <br />
pytz                          2023.3 <br />
PyYAML                        6.0 <br />
pyzmq                         25.0.2 <br />
qtconsole                     5.4.3 <br />
QtPy                          2.3.1 <br />
requests                      2.30.0 <br />
rfc3339-validator             0.1.4 <br />
rfc3986-validator             0.1.1 <br />
scanpy                        1.9.8 <br />
scikit-learn                  1.2.2 <br />
scipy                         1.10.1 <br />
scvelo                        0.3.2 <br />
seaborn                       0.13.2 <br />
Send2Trash                    1.8.2 <br />
session-info                  1.0.0 <br />
setuptools                    67.7.2 <br />
shap                          0.41.0 <br />
sip                           6.7.9 <br />
six                           1.16.0 <br />
sklearn                       0.0.post5 <br />
slicer                        0.0.7 <br />
sniffio                       1.3.0 <br />
soupsieve                     2.3.2.post1 <br />
stack-data                    0.6.2 <br />
statsmodels                   0.14.1 <br />
stdlib-list                   0.10.0 <br />
sympy                         1.12 <br />
terminado                     0.17.1 <br />
texttable                     1.7.0 <br />
threadpoolctl                 3.1.0 <br />
tinycss2                      1.2.1 <br />
toml                          0.10.2 <br />
tomli                         2.0.1 <br />
torch                         2.0.1 <br />
torchaudio                    2.0.2 <br />
torchvision                   0.15.2 <br />
tornado                       6.3 <br />
tqdm                          4.65.0 <br />
traitlets                     5.9.0 <br />
typing_extensions             4.5.0 <br />
tzdata                        2023.3 <br />
umap-learn                    0.5.6 <br />
urllib3                       2.0.2 <br />
wcwidth                       0.2.6 <br />
webencodings                  0.5.1 <br />
websocket-client              1.5.1 <br />
wheel                         0.40.0 <br />
widgetsnbextension            4.0.7 <br />
xgboost                       1.7.6 <br />
zipp                          3.15.0 <br />
