
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


