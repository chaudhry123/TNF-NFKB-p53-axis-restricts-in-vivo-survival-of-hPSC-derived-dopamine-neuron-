
#Add edits for 5sa normalization
rm(list = ls())
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
library("monocle3")
library("org.Hs.eg.db")
addTaskCallback(function(...) {set.seed(100);TRUE})

#Check WT1 and WT2 not switched





#load object Hyunwoo, Needs alignment correcting from your working dir
sce1 <- readRDS("Fayzan.4.SeuratObject.20220420.rds")

fayzan.4 <- sce1
# Take great extra care in counting WT and KO cells!!
# The result should be:
#Always check WT and KO annotations!!
#  fayzan.4@meta.data$Sample %>% table
#WT1  WT2  KO1  KO2
#7342 11235 11730 13129
#fayzan.4@meta.data$condition %>% table
#WT   KO
#18577 24859

# Before correcting,
fayzan.4@meta.data$sample %>% table
#  KO1  KO2  WT1  WT2
# 11235 13129 11730 7342
fayzan.4@meta.data$condition %>% table
#   WT   KO
# 19072 24364

fayzan.4@meta.data$sampleold = fayzan.4@meta.data$sample
fayzan.4@meta.data$sample = "1"
fayzan.4@meta.data$sample[fayzan.4@meta.data$sampleold == "WT2"] = "WT1"
fayzan.4@meta.data$sample[fayzan.4@meta.data$sampleold == "KO1"] = "WT2"
fayzan.4@meta.data$sample[fayzan.4@meta.data$sampleold == "WT1"] = "KO1"
fayzan.4@meta.data$sample[fayzan.4@meta.data$sampleold == "KO2"] = "KO2"
fayzan.4@meta.data$sample = factor(fayzan.4@meta.data$sample, levels=c("WT1", "WT2", "KO1", "KO2"))
fayzan.4@meta.data$sample %>% table
#  WT1  WT2  KO1  KO2
# 7342 11235 11730 13129
fayzan.4@meta.data$condition = factor(gsub("[12]", "", fayzan.4@meta.data$sample), levels=c("WT", "KO"))
fayzan.4@meta.data$condition %>% table
#   WT   KO
# 18577 24859

#Load from your working dir
fayzan.igraph.csv = fread("Fayzan.4.igraph.cutn.7.20220423.csv")
Idents(fayzan.4) = factor(as.character(fayzan.igraph.csv$igraph_cluster), levels=c("1", "2", "3", "4", "5", "6", "7"))

fayzan.4@meta.data$igraph_cluster = factor(as.character(fayzan.igraph.csv$igraph_cluster), levels=c("1", "2", "3", "4", "5", "6", "7"))


combined_sce1 <- as.SingleCellExperiment(fayzan.4)
combined_sce1$clusters <- fayzan.4@meta.data$igraph_cluster
#saveRDS(combined_sce1, "Final.RDS") #Use condition clusters

#Analysis in Seurat
#sampnames <- combined_sce1$condition
#sampnames2 <- substr(sampnames, 0, 2)
combined_sce2 <- combined_sce1

#gene <- rowData(combined_sce2)$Symbol
#temp <- make.unique(gene, sep ="_")
#rownames(combined_sce2) <- temp
sce2 <- as.Seurat(combined_sce2)
#tmpIdent <- sampnames2
Idents(sce2) <- sce2$clusters

#
sce1 <- sce2


sce1
##############
#Monocle3
color_samp <- wes_palette("Darjeeling1",3)
wt_col <- color_samp[2]
ko_col <- color_samp[3]

f1 <- DimPlot(sce1, reduction = "UMAP", group.by = "condition")
print(f1)
expression_matrix <- sce1@assays$originalexp
gene_symbol <- as.list(org.Hs.egSYMBOL)
raw_count_data <- GetAssayData(sce1, assay = "originalexp", slot = "counts")
class(raw_count_data)
cells_info <- sce1@meta.data
gene_name <- rownames(raw_count_data)
gene_name <- sapply(gene_name, function(x) x[[1]][1])

#preparing cds
gene_name <- ifelse(is.na(gene_name), names(gene_name), gene_name) #Check this should be gene name
gene_short_name <- gene_name
gene_id <- rownames(raw_count_data)
genes_info <- cbind(gene_id, gene_short_name)
genes_info <- as.data.frame(genes_info)
rownames(genes_info) <- rownames(raw_count_data)

cds <- new_cell_data_set(expression_data = raw_count_data,
                         cell_metadata = cells_info,
                         gene_metadata = genes_info)

#replace monocle dimensions with seurat
#cds@reducedDims$UMAP <-  sce1@reductions$UMAP@cell.embeddings

cds <- preprocess_cds(cds, num_dim = 50)
#cds <- align_cds(cds, alignment_group = "Sample", alignment_k = 37) #37
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)

#replace monocle clusters with seurat
reducedDims(cds)$UMAP <-  sce1@reductions$UMAP@cell.embeddings
#cds@reducedDims$UMAP <-  sce1@reductions$UMAP@cell.embeddings
cds@clusters$UMAP$partitions <- sce1@meta.data$igraph_cluster
names(cds@clusters$UMAP$partitions) <- rownames(sce1@meta.data)
cds@clusters$UMAP$clusters <- sce1@meta.data$igraph_cluster
names(cds@clusters$UMAP$clusters) <- rownames(sce1@meta.data)


plot_cells(cds)




f2 <- plot_cells(cds,color_cells_by = "partition")
print(f2)
cds <- learn_graph(cds, close_loop=FALSE, use_partition=FALSE)
pdf("monocle_clusters.pdf",width = 6, height = 4)
plot_cells(cds,
           color_cells_by = "igraph_cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE, label_cell_groups =  FALSE, show_trajectory_graph = TRUE)
dev.off()
plot_cells(cds,
           color_cells_by = "Sample",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE, label_cell_groups =  FALSE, show_trajectory_graph = TRUE)

cds <- order_cells(cds)
pdf("monocle_pseudotime.pdf",width = 6, height = 4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()




#Gene Expression Maps
AFD_genes <- c("NR4A2", "TH", "HES5", "MAP2")
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$orig.ident %in% c("igraph_cluster")]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)
pdf("monocle_genes.pdf",width = 6, height = 4)
plot_cells(cds, genes=AFD_genes,
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()
pdf("monocle_genes_baxjun.pdf",width = 6, height = 4)
plot_cells(cds, genes=c(AFD_genes, "BAX", "JUN"),
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()

#Heatmap
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$igraph_cluster)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("igraph_cluster ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=3)
#Save monocle3
save_monocle_objects(cds=cds, directory_path='monocle_surv')

#
sce3 <- subset(sce1, ident==7)
