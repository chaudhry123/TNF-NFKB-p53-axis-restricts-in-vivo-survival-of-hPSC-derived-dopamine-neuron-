
# Bioinformatic methods (NEEDS EDITS)

>TNF NFKB p53 axis restricts in vivo survival of hPSC derived dopamine neuron

Ongoing, first-in-human clinical trials illustrate the feasibility and translational potential of human pluripotent stem cell (hPSC)-based cell therapies in Parkinson’s disease (PD). However, a major unresolved challenge in the field is the extensive cell death following transplantation with <10% of grafted dopamine neurons surviving. Here, we performed a pooled CRISPR/Cas9 screen to enhance survival of postmitotic dopamine neurons in vivo. We identified p53-mediated apoptotic cell death as major contributor to dopamine neuron loss and uncovered a causal link of TNFa-NFκB signaling in limiting cell survival. As a translationally applicable strategy to purify postmitotic dopamine neurons, we performed a cell surface marker screen that enabled purification without the need for genetic reporters. Combining cell sorting with adalimumab pretreatment, a clinically approved and widely used TNFa inhibitor, enabled efficient engraftment of postmitotic dopamine neurons leading to extensive re-innervation and functional recovery in a preclinical PD mouse model. Thus, transient TNFa inhibition presents a clinically relevant strategy to enhance survival and enable engraftment of postmitotic human PSC-derived dopamine neurons in PD.

Here, we provide an overview of the data and scripts that we used for the paper. 

**CONTENT**:

* [Bulk RNA-seq](#Bulk-RNAseq-processingq)
* [CRISPR](#CRISPR-screen-processing)
* [Single Cell](#Single-cell-RNA-seq-Processing)
* [Data for download](#data-for-download)

Scripts were written by Fayzan Chaudhry and Hyein Cho.
Samples were processed at the the sequencing facilities at MSKCC and Weill Cornell Medicine.
Don't hesitate to get in touch with questions related to the code.

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)

## Bulk RNA-seq

### Sample Preparation 



### Processing 




## scVDJ and scRNA-seq data processing

In principle, we followed the ideas and workflows discussed by Amezquita et al. in the excellent online book ["Orchestrating single-cell analysis"](https://osca.bioconductor.org/) describing single-cell analyses carried out with packages of the Bioconductor environment.

### Sample prep

Female NOD mice aged 13-20 weeks were used. Samples were prepared as described above for RNA-seq, except before pooling samples were incubated with hashtags (Biolegend Total-Seq C0301-10) to enable multiplexing.
Samples were stained with NRP-V7 tetramer, Live/dead Zombie dye, and antibodies against CD8alpha and CD45.1; in addition, Biolegend feature barcoding antibodies against CD44, CD62L, CD127, CD73, PD1, CD38, and CD39 were added for CITE-seq analysis. 
2,000-16,000 cells were sorted into PBS+0.04% BSA, immediately analyzed for viability, and processed. 
For details of the library preparation, see [`wetLabPrep_details.md`](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/wetLabPrep_details.md).

### Processing and analysis

VDJ and single-cell RNA-seq including hash-tagged oligos (HTO; used to label cells coming from the same mouse), and antibody-derived tags (ADT; used to label cells based on the expression of selected surface markers) data were processed following the recommendations by [Amezquita et al](https://bioconductor.org/books/release/OSCA) (version: 2020-11-13). In brief, raw read files were processed and aligned using the `CellRanger` pipeline (`cellranger-5.0.0` with `refdata-gex-mm10-2020-A` and `vdj_GRCm38_alts_ensembl-5.0.0`; for annotation supplied by 10X Genomics: <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>).
Read counts representing gene expression were used for removing barcodes representing empty droplets or droplets with high amounts of mitochondrially encoded gene products: using the `isOutlier()` function of the scater package, all cells whose fractions of mitochondrial reads were higher than 3x median absolute deviation (MAD) were flagged and removed as well as cells with fewer than 10^2.5 detected genes.
Following the first round of filtering, HTO demultiplexing and doublet detection were done with [`CiteFuse`](http://www.bioconductor.org/packages/release/bioc/html/CiteFuse.html), removing more droplets and identifying the donor mouse for each remaining droplet. After removal of low-quality droplets and barcodes associated with HTO-inferred doublets, ADT values were normalized using
`DropletUtils::inferAmbience()` and `scuttle::medianSizeFactors()`.
Finally, droplets for which less than 50% of the ADT were detected compared to other droplets were removed and only genes that were expressed in at least 5 cells per sample were kept for downstream analyses.

VDJ sequencing data was imported into R using the `import_vdj()` function of the [`djvdj` package](https://github.com/rnabioco/djvdj), which was adjusted to work with `SingleCellExperiment` objects (<https://github.com/friedue/SCEdjvdj>).

Size-factor normalized logcounts were obtained via `scran::computeSumFactors()` and `scater::logNormCounts()` (Lun2016, McCarthy2017).
Cells from the three different technical replicates were integrated with `batchelor::fastMNN()` using the top 2500 most variable genes with min. mean normalized expression of 0.001 (Haghverdi2018).
Clustering was performed with `igraph::cluster_louvain()`, dimensionality reductions were done with `scater` functions (`runUMAP()`) using the batch-corrected values and `destiny::diffusionMap()` (Lun2016clustering, Angerer2015).
Marker genes were determined with `scran::findMarkers()`.
The `TSCAN package (v.1.28.0)` was used to calculate pseudotime values and trajectories as well as genes associated with the pseudotime gradients.
GO term enrichments were calculated with the `clusterProfiler` package's functions `compareCluster()` and `enrichGO()` after excluding ribosomal genes from the gene lists of interest.
All plots were generated using `ggplot2` packages and the `pheatmap` package for heatmaps.

For details, see [processing_scRNAseq.md](scRNAseq/processing_scRNAseq.md), [processingVDJseq.md](scRNAseq/processingVDJseq.md), [figures_scRNAseq.Rmd](scRNAseq/figures_scRNAseq.Rmd) and [figures_VDJseq.Rmd](scRNAseeq/figures_VDJseq.Rmd).

## Integration with public scRNA-seq data sets

Gene count matrices for day 7 CD8 T cells from acute and chronic infection were obtained from GEO (GSE119940); using cell labels provided by Chen Yao we extracted the data for barcodes corresponding to memory precursor and memory-like cells as described in Yao et al. (2019).
scRNA-seq data for Tcm were downloaded from GEO (day 129 following acute infection: GSM3732587) where they had been uploaded to by Schauder et al.
scRNA-seq of Tcms, MPECS, and Tpex were subsequently integrated with the pLN data set using `batchelor::multiBatchNorm()` and `batchelor::fastMNN()` as described above.
For global comparisons of the different populations of T cells, we created pseudo-bulk samples by aggregating the read counts per gene across cells of the same population.
These were then cpm-normalized via `edgeR::calcNormFactors` (Robinson 2010) and subsequently analyzed and visualized via PCA and hierarchical clustering using base R functions as well as `pcaExplorer::hi_loadings()` and the `dendextend` package (Marini 2019, Galili 2015).
All other analyses were done with the same principles and packages as described above.

For details see [Schauder2021.Rmd](scRNAseq/Schauder2021.Rmd), [Yao2019.Rmd](scRNAseq/Yao2019.Rmd), [figures_public_scRNAseq.Rmd](scRNAseq/figures_public_scRNAseq.Rmd) and the corresponding PDF/HTML files.



## Package versions

## Data for download

The raw data (fastq files, read counts from CellRanger) can be downloaded from GEO (INSERT NUMBER HERE).


>The easiest way to get started is to use the processed data provided here.

For the single-cell data, some of the data can be downloaded from Box in the form of RDS (load into R via `in_data <- readRDS()`) or RDA objects (load into R via `load()`).
For the direct links to the RDS/RDA objects, see [`scRNAseq/data/data_links.txt`](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/data/data_links.txt).
The `data/` directory in the scRNA-seq directory contains some text files that contain just the cell labels and the mouse labels for individual cells.


