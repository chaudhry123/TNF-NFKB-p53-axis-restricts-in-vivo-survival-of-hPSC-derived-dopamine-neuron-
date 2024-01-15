
# Bioinformatic methods

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




## scRNA-seq data processing

In principle, we followed the ideas and workflows discussed by Amezquita et al. in the excellent online book ["Orchestrating single-cell analysis"](https://osca.bioconductor.org/) describing single-cell analyses carried out with packages of the Bioconductor environment.

### Sample prep


### Processing and analysis





## Package versions

## Data for download

The raw data (fastq files, read counts from CellRanger) can be downloaded from GEO (INSERT NUMBER HERE).


>The easiest way to get started is to use the processed data provided here.

For the single-cell data, some of the data can be downloaded from Box in the form of RDS (load into R via `in_data <- readRDS()`) or RDA objects (load into R via `load()`).
For the direct links to the RDS/RDA objects, see [`scRNAseq/data/data_links.txt`](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/data/data_links.txt).
The `data/` directory in the scRNA-seq directory contains some text files that contain just the cell labels and the mouse labels for individual cells.


