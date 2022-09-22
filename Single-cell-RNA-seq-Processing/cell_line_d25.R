#Import libraries
rm(list = ls())
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
library(Seurat)
library(BiocParallel)
library(Matrix)

addTaskCallback(function(...) {set.seed(100);TRUE}) #Permanently sets seed so is reproducible

#Load in data of 3 wildtypes and 2 knockouts
#path1="WT1outs/filtered_feature_bc_matrix" 
#path2="WT2outs/filtered_feature_bc_matrix"
#path3="WT3outs/filtered_feature_bc_matrix"
#path4="KO1outs/filtered_feature_bc_matrix"
#path5="KO2outs/filtered_feature_bc_matrix"
path6 = "/Users/fayzan/outs/filtered_feature_bc_matrix"
combined_sce1 <- read10xCounts(c(path6), c("cellline"))
#combined_sce2 <- read10xCounts(path1, "WT1")
#combined_sce1 <- read10xCounts(path1, "WT1")
#samp <- cbind(combined_sce2, combined_sce)

#Map Mito
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(combined_sce1)$ID, column="SEQNAME", keytype="GENEID")
#combined_sce3 <- combined_sce
#QC
min_UMI <- 1000
min_nGene <- 500
max_nGene <- 7000
max_mito_frac <- 0.25

combined_sce1 <- addPerCellQC(combined_sce1, subsets=list(Mt=which(location=="MT")), flatten=TRUE)
combined_sce1 <- addPerFeatureQC(combined_sce1, flatten=TRUE)
#samples <- unique(combined_sce1$Sample)

##plot number of UMI
jpeg("jpeg1.jpeg")
plotColData(combined_sce1, x="Sample", y="sum", colour_by="Sample") +
  ylab("Number of UMIs") + geom_hline(yintercept=min_UMI, color='red')
dev.off()

## plot number of genes
jpeg("jpeg2.jpeg")
plotColData(combined_sce1, x="Sample", y="detected", colour_by="Sample") +
  ylab("Number of genes")+ geom_hline(yintercept=min_nGene, color='red') + geom_hline(yintercept=max_nGene, color='red')
dev.off()

## plot number of Mt. counts
jpeg("jpeg3.jpeg")
plotColData(combined_sce1, x="Sample", y="subsets_Mt_sum", colour_by="Sample") + ylab("Number of Mt UMIs")
dev.off()

## plot % of Mt. counts
jpeg("jpeg4.jpeg")
plotColData(combined_sce1, x="Sample", y="subsets_Mt_percent", colour_by="Sample") + ylab("% of Mt UMIs") +
  geom_hline(yintercept=max_mito_frac*100, color='red')
dev.off()


## histogram of total UMI counts
jpeg("jpeg5.jpeg")
ggplot(data=tibble::as_tibble(colData(combined_sce1)), aes(x=sum)) + geom_histogram(binwidth=100, color='blue') +
  facet_wrap(~ Sample, ncol=5) + xlab("Number of UMIs") + geom_vline(xintercept=min_UMI, color='red')
dev.off()

## histogram of total gene counts
jpeg("jpeg6.jpeg")
ggplot(data=tibble::as_tibble(colData(combined_sce1)), aes(x=detected)) + geom_histogram(binwidth=50, color='lightblue') +
  facet_wrap(~ Sample, ncol=5) + xlab("Number of genes") + geom_vline(xintercept=min_nGene, color='red')
dev.off()

## number of genes vs. number of UMI color by %Mt genes
jpeg("jpeg7.jpeg")
ggplot(data=tibble::as_tibble(colData(combined_sce1)), aes(x=detected, y=sum)) + geom_point(aes(col=subsets_Mt_percent)) +
  facet_wrap( ~Sample, ncol=5) + xlab("Gene counts") + ylab("UMI counts")
dev.off()

jpeg("jpeg8.jpeg")
count.depth <- tibble::as_tibble(colData(combined_sce1)) %>% arrange(desc(sum))
ggplot(data= count.depth, aes(x=1:nrow(count.depth), y=total)) + geom_point(aes(col=total)) + facet_wrap( ~Sample, ncol=5)
dev.off()

jpeg("jpeg13.jpeg")
count.depth <- tibble::as_tibble(colData(combined_sce1)) %>% arrange(desc(sum))
ggplot(data= count.depth, aes(x=1:nrow(count.depth), y=total)) + geom_point(aes(col=total)) #+ facet_wrap( ~Sample, ncol=2)
dev.off()

## filter by UMI and gene counts
print(paste("Size of before cell filtering", ncol(combined_sce1), "cells", nrow(combined_sce1), "genes"))
out.of.range.UMI <- (combined_sce1$total <= min_UMI ) 
out.of.range.genecounts <- (combined_sce1$detected <= min_nGene |  combined_sce1$detected >= max_nGene)

combined_sce1 <- combined_sce1[,(!out.of.range.UMI & !out.of.range.genecounts)]
print(paste("Size of matrix after cell filtering", ncol(combined_sce1), "cells", nrow(combined_sce1), "genes"))

## Filter by %Mt
high.mt <- (colData(combined_sce1)$subsets_Mt_percent >= (max_mito_frac*100))
combined_sce1 <- combined_sce1[,!high.mt]
print(paste("Size of matrix after Mitochondria filtering", ncol(combined_sce1), "cells", nrow(combined_sce1), "genes"))


## remove low expressed genes, minimumly expressed 
gene2keep <- rowData(combined_sce1)$detected > 0

#PLOT
jpeg("jpeg14.jpeg")
plot(scater::nexprs(combined_sce1, byrow=TRUE), rowData(combined_sce1)$detected, pch=19, col=ifelse(gene2keep, "gray60", "red"),
     xlab="number of cells with non-zero counts", ylab="% dropout")
abline(h=quantile(rowData(combined_sce1)$detected, prob=0.01), col='red')
dev.off()

combined_sce1 <- combined_sce1[gene2keep,]
print(paste("Size of matrix after gene dropout removal", ncol(combined_sce1), "cells", nrow(combined_sce1), "genes")) 

#PLOT
jpeg("jpeg15.jpeg")
scater::plotRowData(combined_sce1, x = "detected", y = "mean") + xlab("percentage of cells with counts for each gene") +
  ylab("mean count of gene") + ggtitle("% of cells gene is detected")
dev.off()

## normalize
clusters <- scran::quickCluster(combined_sce1, min.size=100)
combined_sce1 <- scran::computeSumFactors(combined_sce1, cluster=clusters)
combined_sce1 <- scater::logNormCounts(combined_sce1)

BinomialDeviance <- function(x, p, n)
  
{
  ## x - gene count vector
  ## p - colSums(x)
  ## n - colSums(x).
  term1 <- sum(x*log(x/(n*p)), na.rm=TRUE)
  nx <- n-x
  term2 <- sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
  return(2*(term1+term2))
  
}




combined_sce1 <- combined_sce1
var.genes <- modelGeneVar(combined_sce1) #STOPPED HERE
#jpeg("jpeg16.jpeg")
#plot(var.genes$mean, var.genes$total, pch=19, col="gray", xlab="Mean lob-expression", ylab="Variance")
#curve(metadata(var.genes)$trend(x), col="blue", add=TRUE)
#dev.off()
#Combine
#combined_sce1 


cnts <- counts(combined_sce1)
cnts <- cnts[Matrix::rowSums(cnts) >0,]
size.fac <- Matrix::colSums(cnts)
BinomialDeviance(cnts[34,], sum(cnts[34,])/sum(size.fac),size.fac)


top.hvgs <- getTopHVGs(var.genes, n=4000)
combined_sce1 <- scater::runPCA(combined_sce1, subset_row = top.hvgs)
#percent.var <- attr(reducedDim(combined_sce1), "percentVar")
#chosen.elbow <- PCAtools::findElbowPoint(percent.var)
#chosen.elbow
#percent.var
combined_sce1 <- denoisePCA(combined_sce1, var.genes)

#numReducedDim <- maxLikGlobalDimEst(as.matrix(logcounts(combined_sce1)), k=10)
combined_sce1 <- scater::runPCA(combined_sce1,  ncomponents = 100, ntop = 4000)

#Plot PCA
jpeg("cellline_pca.jpeg")
plotReducedDim(combined_sce1, dimred="PCA", ncomponents=4,
               colour_by="Sample") + geom_point(aes(size = .1, colour = combined_sce1$Sample))
dev.off()

#Correct for batch effects
#merged <- correctExperiments(combined_sce1, 
#                             batch=combined_sce1$Sample, 
#                             subset.row=top.hvgs,
#                             PARAM=FastMnnParam(
#                               merge.order=list(
#                                 list(1,2,3), # WT (3 replicates)
#                                 list(4,5)  # KO (2 replicates)
#                               )
#                             )
#)

## clustering
g <- buildSNNGraph(combined_sce1, k=40, use.dimred="PCA")
# clusters <- igraph::cluster_walktrap(g)$membership
temp2 <- igraph::cluster_walktrap(g)
#saveRDS(combined_sce1, file="umappca100.rds")
#saveRDS(temp2, "temp2.rds")
#saveRDS(g, "g.rds")
#combined_sce1 <- readRDS("final_umap.rds")
temp1 <- igraph::cut_at(temp2, n=3)
ratio <- pairwiseModularity(g, temp1, as.ratio=TRUE)
#dim(ratio)
#jpeg("cluster2.jpeg")
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE, color=colorRampPalette(c("white", "blue"))(100))
#dev.off()
#saveRDS(temp, "temp_cell.rds")
#saveRDS(g, "g_cell.rds")
#saveRDS(combined_sce1, "umappca100_temp.rds")
#g <- readRDS("g_cell.rds")
#combined_sce1  <- readRDS("umap_cellline.rds")
#temp2 <- readRDS("temp_cell.rds")
combined_sce1$cluster <-  factor(temp1)
sce1 <- readRDS("celllinevsko.rds")
combined_sce1 <- scater::runUMAP(combined_sce1, dimred="PCA", min_dist=.01, n_neighbors = 50)
#saveRDS(combined_sce1, file="umap_cellline.rds")

#combined_sce1 <- readRDS("new_umap_4samp.rds")
jpeg("jpeg17_1.jpeg")
scater::plotUMAP(combined_sce1, colour_by="cluster", text_by="cluster")
dev.off()
#saveRDS(combined_sce1, "finalumap2.rds")
#KEeP at 7
#CLuster metrics
#sil.approx <- approxSilhouette(reducedDim(combined_sce1, "PCA"), clusters=colLabels(combined_sce1))
#sil.data <- as.data.frame(sil.approx)
#sil.data$closest <- factor(ifelse(sil.data$width > 0, colLabels(combined_sce1), sil.data$other))
#sil.data$cluster <- colLabels(combined_sce1)

#ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
 # ggbeeswarm::geom_quasirandom(method="smiley")

#PURITY
#pure.pbmc <- neighborPurity(reducedDim(combined_sce1, "PCA"), colLabels(combined_sce1))
#pure.pbmc
#pure.data <- as.data.frame(combined_sce1)
#pure.data$maximum <- factor(pure.data$maximum)
#pure.data$cluster <- colLabels(combined_sce1)

#ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
#  ggbeeswarm::geom_quasirandom(method="smiley")

#GRAPH
cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1),
                                                  mode="upper", weighted=TRUE, diag=FALSE)

# Increasing the weight to increase the visibility of the lines.
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*1,
     layout=igraph::layout_with_lgl)
