## Single-cell-RNA-seq-Processing

### Description

To define the transcriptional landscape of grafted dopamine neurons, we performed single cell mRNA sequencing from p53 wild-type (WT) and p53 knock-out (KO) grafted neurons, re-isolated from the mouse brain at 1 dpt. This is the processing pipeline to generate the figures used in the accompanying paper. For full transparency and reproducibility, we've provided the code implementation, including scripts, data, and dependencies.

* This repository includes a comprehensive breakdown of each step in the processing pipeline, from data preprocessing to figure generation. We've organized the code into logical modules to facilitate understanding and modification. Detailed comments and documentation are included to guide you through the implementation and parameter tuning process.

* To ensure that you can replicate our results precisely, we've also included the datasets used in our experiments. These datasets are properly formatted and annotated, making them ready for direct use with the provided scripts. The R data objects for this repository are available to download at our Box foler (https://mskcc.box.com/s/wn5uvwxu2xm4hw219mo0id3r9nkyyprx).

* If you would like to generate these `.rds` files independently we included the preprocessing steps applied to the raw data, enabling you to reproduce our preprocessing pipeline if needed.

* The codebase is designed to be modular and extensible, allowing for easy adaptation to new datasets or modifications to existing processing steps.

* We included instructions on how to set up the necessary environment to run the code, including information on required dependencies and how to install them. We recommend using virtual environments to isolate the dependencies for this project, ensuring compatibility and reproducibility across different systems.

* Please [email](ffc4001@med.cornell.edu) us for bug reports or with suggestion for optimization. Feel free to open an issue or submit a pull request on GitHub.

### Executing program
* Pre-processed single cell data objects can be downloaded from [this Box folder](https://mskcc.box.com/s/wn5uvwxu2xm4hw219mo0id3r9nkyyprx), or they can be generated using the scripts provided here. If you choose to use the downloaded objects you may need to modify the code where these files are loaded (i.e. `readRDS` calls) to match your own working directory. For example, in `Figure1.spin.Rmd` the files are assumed to be in `./TNF-NFKB-2024/` relative to the current working directory. 

* To obtain the processed single-cell experiment from the combined day 1 replicates, execute the `combinedday1replicates.Rmd` R markdown. Subsequently, apply additional processing steps by running `processing.R` that generates the `FinalObj.rds` file. 

* For day 25 in-vitro cell line dataset, run the `cell_line_d25.Rmd` script to generate the processed single-cell experiment (`cellline_clean.rds`) file.

* After obtaining the necessary objects, run the `Figures1.spin.Rmd` script to convert them into a `Seurat` objects and generate the figures for Figure 5 and Supplement Figure 5. The required objects for this step are the `cellline_clean.rds` and `FinalObj.RDS`.

* For Monocle trajectory analysis, execute the `Monocle.R` script and the `Jupyter` Python notebook for survival analysis that generates the Monocle results for the supplement.

### Authors
Fayzan Chaudhry ffc4001@med.cornell.edu
