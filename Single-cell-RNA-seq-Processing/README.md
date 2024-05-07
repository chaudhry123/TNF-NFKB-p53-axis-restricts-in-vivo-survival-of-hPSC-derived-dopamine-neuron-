# Single-cell-RNA-seq-Processing

## Description

To further define the transcriptional landscape of grafted dopamine neurons, we performed single cell mRNA sequencing from p53 wild-type (WT) and p53 knock-out (KO) grafted neurons, re-isolated from the mouse brain at 1 dpt. This is the processing pipeline to generate the figures used in the accompanying paper.For full transparency and reproducibility, we've provided the code implementation, including scripts, data, and dependencies.

In this repository, you'll find a comprehensive breakdown of each step in the processing pipeline, from data preprocessing to figure generation. We've organized the code into logical modules to facilitate understanding and modification. Detailed comments and documentation are included to guide you through the implementation and parameter tuning process.

To ensure that you can replicate our results precisely, we've also included the datasets used in our experiments. These datasets are properly formatted and annotated, making them ready for direct use with the provided scripts. Additionally, we've documented any preprocessing steps applied to the raw data, enabling you to reproduce our preprocessing pipeline if needed.

The codebase is designed to be modular and extensible, allowing for easy adaptation to new datasets or modifications to existing processing steps. Whether you're interested in reproducing our results, applying our methods to your own data, or extending our work with new techniques, the repository provides a solid foundation to build upon.

Furthermore, we've included instructions on how to set up the necessary environment to run the code, including information on required dependencies and how to install them. We recommend using virtual environments to isolate the dependencies for this project, ensuring compatibility and reproducibility across different systems.

In addition to the code and data, we've provided documentation on the experimental setup, including parameter configurations and hardware specifications used during testing. This information is crucial for understanding the context in which the figures were generated and for reproducing our experiments under similar conditions.

We encourage collaboration and welcome contributions from the community to further improve the codebase and expand its functionality. Whether you spot a bug, have a suggestion for optimization, or want to propose a new feature, feel free to open an issue or submit a pull request on GitHub.

By sharing our code and data openly, we aim to foster transparency, reproducibility, and collaboration in scientific research, ultimately advancing our collective understanding of the subject matter. We hope you find this resource useful and look forward to any feedback or contributions you may have.

### Executing program

To obtain the processed single-cell experiment from the combined day 1 replicates, execute the combinedday1replicates script. Subsequently, apply additional processing steps using the provided processing script. For the day 25 in vitro cell line dataset, run the cell_line_d25 script to generate the processed single-cell experiment.

After obtaining the necessary objects, run the Figures1 script to convert them into a Seurat object and generate the figures for Figure 5 and Supplement Figure 5. The required objects for this step are the cell line object and finalobj.rds.

For Monocle analysis, execute the Monocle.R script and the Jupyter notebook for survival analysis to generate Monocle results for the supplement.

These steps ensure the generation of processed single-cell experiments, conversion to Seurat objects, and the creation of figures for analysis and visualization, along with additional analysis using Monocle and survival analysis.

## Authors

Contributors names and contact info

Fayzan Chaudhry Fachaudhry96@gmail.com
