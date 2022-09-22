library("SingleR")
library("celldex")
library("scRNAseq")
combined_sce1 <- readRDS("umappca100.rds")
ref.data1 <- scRNAseq::LaMannoBrainData("mouse-adult")
#libsizes <- colSums(counts)
#size.factors <- libsizes/mean(libsizes)
#logcounts(ref.data1) <- log2(t(t(counts)/size.factors) + 1)
hpca.se <- HumanPrimaryCellAtlasData()
pred.hesc <- SingleR(test = sceM1, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
predictions <- SingleR(test=combined_sce1, assay.type.test=2, 
                       ref=sce, labels =sce$`level1 class`, clusters=combined_sce1$cluster)
plotScoreHeatmap(predictions, clusters=unique(ref.data1))

combined_sce1 <- scater::runUMAP(combined_sce1, dimred="PCA", min_dist=.01, n_neighbors = 50)


scater::plotUMAP(combined_sce1, colour_by="cell", text_by="cell")

current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("Mature Dopamine A9", "Mature Dopamine", "Mature Dopamine", "Progenitor", "Low in Map2/Unclear", "Mature Dopamine A9 and A10", "Progenitor")
combined_sce1$cell <- plyr::mapvalues(x = combined_sce1$cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "tsne", do.label = TRUE, pt.size = 0.5)
scater::plotUMAP(combined_sce1, colour_by="cell", text_by="cell", cols ="red")

#LEMANNO ANNOTATION
sce2 <- scRNAseq::LaMannoBrainData("human-embryo")


# SingleR() expects reference datasets to be normalized and log-transformed.
library(scuttle)

sce2 <- logNormCounts(sce2) 
#sceG <- sceG[,1:100]




sce
#de_genes <- unique(c(mDA.mature.markers, mDA.precursor.markers, A10.markers, A9.markers, non.mDA.markers))
predictions2 <- SingleR(test=combined_sce1, assay.type.test=2, 
                       ref=sce3, labels =sce3$Cell_type, de.method="wilcox")
png("Feb_Plots/WT2_assignment_final_final3.png")
plotScoreHeatmap(predictions2, clusters=unique(combined_sce1$clusters))
dev.off()
current.cluster.ids <- c("1", "2")
new.cluster.ids <- c("hNbM", "hProgFPL")
combined_sce1$cell <- plyr::mapvalues(x = combined_sce1$cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "tsne", do.label = TRUE, pt.size = 0.5)
scater::plotUMAP(combined_sce1, colour_by="cell", text_by="cell")

#Find marker genes
by.t <- scran::pairwiseTTests(assay(sce, 2), sce$Cell_type, direction="up")
markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=50)
trained <- trainSingleR(sce, sce$Cell_type, genes=markers)
length(trained$common.genes)
predictions3 <- SingleR(test=combined_sce1, assay.type.test=1, 
                        ref=sce, labels=sce$Cell_type, genes=markers)
current.cluster.ids <- c("1", "2")
new.cluster.ids <- c("hNbML1", "hProgFPL")
combined_sce1$Cell <- predictions2$pruned.labels
DimPlot(object = pbmc, reduction = "tsne", do.label = TRUE, pt.size = 0.5)
combined_sce1$cluster1 <- predictions8$pruned.labels
scater::plotUMAP(combined_sce1, colour_by="Cell")  #+ scale_color_manual(values=mycolors)

#HERE
library(scuttle)
#sce <- scRNAseq::LaMannoBrainData("human-ips", ensembl = TRUE)
sce <- logNormCounts(sce) 
predictions9 <- SingleR(test=combined_sce1,
                        ref=sce, labels=sce$Cell_type, de.method="wilcox",assay.type.test=2, clusters=combined_sce1$clusters)
plotScoreHeatmap(predictions9$pruned.labels)
table <- data.frame(unclass(table(predictions7$pruned.labels)))
table$unclass.table.predictions7.pruned.labels.. <- as.numeric(table$unclass.table.predictions7.pruned.labels..)
pie(labels = rownames(table), x =table$unclass.table.predictions7.pruned.labels..,
    main="Pie Chart of Cell Annotations")



#Pie chart
c25 <- c(
  "dodgerblue2", "#E31A1C", # red

  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "yellow4", "yellow3",
  "darkorange4", "brown"
)
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(22)
ggplot(data = table, aes(x = "", y = table$unclass.table.predictions7.pruned.labels.., fill = rownames(table))) + 
  geom_bar(stat = "identity", color = "black") +
  labs(fill = "Cell Type") +
  coord_polar("y") +
  theme_void() +scale_fill_manual(values = mycolors)

#load object Hyunwoo
sce1 <- readRDS("Fayzan.4.SeuratObject.20220420.rds")
fayzan.4 <- sce1
# Take great extra care in counting WT and KO cells!!
# The result should be:
Always check WT and KO annotations!!
fayzan.4@meta.data$Sample %>% table
WT1  WT2  KO1  KO2
7342 11235 11730 13129
fayzan.4@meta.data$condition %>% table
WT   KO
18577 24859

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
library(data.table)
fayzan.igraph.csv = fread("Fayzan.4.igraph.cutn.7.20220423.csv")
Idents(fayzan.4) = factor(as.character(fayzan.igraph.csv$igraph_cluster), levels=c("1", "2", "3", "4", "5", "6", "7"))

fayzan.4@meta.data$igraph_cluster = factor(as.character(fayzan.igraph.csv$igraph_cluster), levels=c("1", "2", "3", "4", "5", "6", "7"))


combined_sce1 <- as.SingleCellExperiment(fayzan.4)
combined_sce1$clusters <- fayzan.4@meta.data$igraph_cluster
sce2 <- scRNAseq::LaMannoBrainData("human-embryo")


# SingleR() expects reference datasets to be normalized and log-transformed.
library(scuttle)

sce2 <- logNormCounts(sce2) 

#de_genes <- unique(c(mDA.mature.markers, mDA.precursor.markers, A10.markers, A9.markers, non.mDA.markers))
sce3 <- subset(sce2 , ,Cell_type == c('hNbM') | Cell_type == c('hNbML1') | Cell_type == c('hProgFPL'))
predictions2 <- SingleR(test=combined_sce1, assay.type.test=2, 
                        ref=sce3, labels =sce3$Cell_type, de.method="wilcox")
combined_sce1$pruned_labels <- predictions2$pruned.labels

sce1 <- as.Seurat(combined_sce1)

DimPlot(sce1)
#LabelClusters(plot, predictions9$pruned.labels)
fayzan.4 <- RenameIdents(object = fayzan.4, `1` = "hNbM", `2` = "hNbM", `3` = "hProgFPL", '4'= "hNbM", "5" = "hProgFPL", "6" = "hProgFPL", "7" = "hPeric")


Idents(sce) <- sce$Cell_type
#sce3 <- subset(x = sce, idents = c('hNbM','hNbML', 'hProgFPL'))
sce3 <- subset(sce2 , ,Cell_type == c('hNbM') | Cell_type == c('hNbML1') | Cell_type == c('hProgFPL'))
sce7 <- subset(sce1 , condition == "KO")
Idents(sce7) <- sce7$pruned_labels
EnhancedVolcano(markers_ko,
                lab = rownames(markers_ko),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = genesofinterest,
                title = 'Top DE Genes hNbM vs hNbML1 KO Only',
                pCutoff = .01,
                FCcutoff = .10,
                pointSize = 3.0,
                labSize = 4.0,
                col=c('blue', 'green', 'black', 'red3'),
                colAlpha = 1,
                xlim =c(-1,1))


markers_ko <- FindMarkers(sce7, logfc.threshold= .05, min.pct = 0.05, min.diff.pct = 0.02, ident.1 ="hNbM", ident.2="hNbML1", test.use = "MAST")
markers_ko <- markers_ko %>% dplyr::filter(avg_log2FC > .05) %>% dplyr::filter(p_val < .01)









#ACE
sce3 <- subset(fayzan.4,  condition == "KO")
DimPlot(sce3, group.by="igraph_cluster")
FeaturePlot(fayzan.4, feature="HES5", order=T, pt.size=3.0) + labs(title="HES5") + scale_color_gradientn(colors=c("grey95", rev(magma(99))), limits=c(0, 5)) 
