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
#load object Hyunwoo, Needs alignment correcting
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

fayzan.igraph.csv = fread("Fayzan.4.igraph.cutn.7.20220423.csv")
Idents(fayzan.4) = factor(as.character(fayzan.igraph.csv$igraph_cluster), levels=c("1", "2", "3", "4", "5", "6", "7"))

fayzan.4@meta.data$igraph_cluster = factor(as.character(fayzan.igraph.csv$igraph_cluster), levels=c("1", "2", "3", "4", "5", "6", "7"))


combined_sce1 <- as.SingleCellExperiment(fayzan.4)
combined_sce1$clusters <- fayzan.4@meta.data$igraph_cluster
saveRDS(combined_sce1, "FinalObj.RDS") #Use condition clusters
