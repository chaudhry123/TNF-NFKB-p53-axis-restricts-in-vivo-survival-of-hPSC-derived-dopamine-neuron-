# Bulk RNA-seq Processing

### Executing program

To obtain the processed bulk RNA-seq experiment, refer to this repository or the GEO repository.

Subsequently, run the file bulk-rnaseq-run-through.sh. It internally invokes R code in bulk-rnaseq-run-through.R.

The code downloads the data and generates the PCA plot, the MA plot, the heat map, and functional enrichment plots.

These steps ensure the generation of processed single-cell experiments, conversion to Seurat objects, and the creation of figures for analysis and visualization, along with additional analysis using Monocle and survival analysis.

## Authors

Contributors names and contact info

Hyein Cho choh@mskcc.org
