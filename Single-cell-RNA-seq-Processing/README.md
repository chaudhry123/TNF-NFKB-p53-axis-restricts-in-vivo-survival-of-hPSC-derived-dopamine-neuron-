# Project Title

Single-cell-RNA-seq-Processing

## Description

To further define the transcriptional landscape of grafted dopamine neurons, we performed single cell mRNA sequencing from p53 wild-type (WT) and p53 knock-out (KO) grafted neurons, re-isolated from the mouse brain at 1 dpt. To further define the transcriptional landscape of grafted dopamine neurons, we performed single cell mRNA sequencing from p53 wild-type (WT) and p53 knock-out (KO) grafted neurons, re-isolated from the mouse brain at 1 dpt. This is the processing pipeline to generate the figures used in the accompanying paper.


### Dependencies
library(DropletUtils)
library(scater)
library(scran)
library(ggplot2)
library(Seurat)
library(dplyr)
library(viridis)
library(EnhancedVolcano)
library(formattable)
library(pheatmap)
library(ape)
library(ComplexHeatmap)
library(patchwork)
library(scales)
library(wesanderson)
library(ggrepel)
library(ggeasy)
library(data.table)
library(scuttle)
library(SingleR)
library(DropletUtils)
library(scater)
library(AnnotationHub)
library(scran)
library(DropletTestFiles)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(ggplot2)
library(Seurat)
library(intrinsicDimension)
library(batchelor)
library(dplyr)
library(bluster)
library(pheatmap)
library(rtracklayer)
library(stringr)
library(data.table)
library(Rsubread)


### Executing program

Run combinedday1replicates.R to obtain the proccessed single cell experiment. Also, run cell_line_d25.R to obtain a processed single cell experiment for day 25 in vitro cellline dataset. Then run Figures.R to convert to a seurat object and obtain the figures.

## Authors

Contributors names and contact info

Fayzan Chaudhry Fachaudhry96@gmail.com
