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
addTaskCallback(function(...) {set.seed(100);TRUE})
#Check WT1 and WT2 not switched





#load object Hyunwoo
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
#saveRDS(combined_sce1, "Final.RDS") #Use condition clusters

#Analysis in Seurat
#sampnames <- combined_sce1$condition
#sampnames2 <- substr(sampnames, 0, 2)
combined_sce2 <- combined_sce1

gene <- rowData(combined_sce2)$Symbol
temp <- make.unique(gene, sep ="_")
rownames(combined_sce2) <- temp
sce2 <- as.Seurat(combined_sce2)
#tmpIdent <- sampnames2
Idents(sce2) <- sce2$clusters

#
sce1 <- sce2


sce1
##############
#UMAP
color_samp <- wes_palette("Darjeeling1",3)
wt_col <- color_samp[2]
ko_col <- color_samp[3]
pdf("Figure4a.pdf", width = 4, height = 4)
my_color_palette <- hue_pal()(length(levels(sce1$clusters)))
DimPlot(sce1)
dev.off()
pdf("Figure4b.pdf", width = 4, height = 4)
DimPlot(sce1, group.by = "condition", cols = c(wt_col, ko_col), label = FALSE) + ggtitle("") + theme(legend.text=element_text(size=13))
dev.off()

#MAKE HEAT MAP, CLUSTERPLOTS, VLN
markers <- FindAllMarkers(sce1, min.pct = 0.15, min.diff.pct = 0.04, logfc.threshold = .25, test.use = "MAST") #Try with MAST
#top20 <- markers %>% dplyr::group_by(markers$cluster) %>% top_n(5, avg_log2FC)
#temp3 <- as.matrix(ScaleData(object = mat, block.size=100))
#png("heatmaptopgroup1month_cluster.png", width = 800, height = 800)
top20 <- markers %>%  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#sce1$
sce1 <- ScaleData(object = sce1, block.size=100)
#png("heatmaptopgroup.png", width = 800, height = 800)
temp3 <- GetAssayData(sce1, slot = "scale.data")
mat<- temp3[unique(c("JUN", top20$gene)), ] %>% as.matrix()
#mat<- t(scale(t(mat)))
heatmapanno <- HeatmapAnnotation(Group = sce1$condition, Cluster = Idents(sce1), col = list(Group=c("WT" = wt_col, "KO" = ko_col), Cluster =c("1" = "#F8766D", "2" =  "#C49A00", "3" =  "#53B400", "4" =  "#00C094", "5" =  "#00B6EB", "6" =  "#A58AFF", "7" ="#FB61D7")))
col_fun = circlize::colorRamp2(c(-1, 0, 2), c("#140b34" , "#84206b", "#f6d746"))
#Figure 4D 
pdf("Figure4d_ALTMAST.pdf", width = 6, height = 6.5)
Heatmap(mat, name = "Scaled Expression",  
        column_split = Idents(sce1),
        #column_title = "Cluster",
        cluster_columns = FALSE, 
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        column_title_gp = gpar(fontsize = 7),
        column_gap = unit(.6, "mm"),
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        row_dend_width = unit(20, "mm"),
        col = col_fun,
        row_names_gp = gpar(fontsize = 8),
        column_title_rot = 0,
        top_annotation = heatmapanno,
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)
dev.off()

#Figure 4F #Exclude 7
sce3 <- subset(sce1, ident = 7, invert = TRUE)
pdf("Figure4F_2.pdf", width = 3, height = 6)
v1 <- VlnPlot(sce3, "BAX", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col), pt.size = 0)
v2 <- VlnPlot(sce3, "TNFRSF12A", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v3 <- VlnPlot(sce3, "JUN", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v1+v2 + v3
dev.off()

pdf("Figure4F_ALT_2.pdf", width = 3, height = 6)
v1 <- VlnPlot(sce3, "BAX", split.by = "condition", cols = c(wt_col, ko_col), pt.size = 0)
v2 <- VlnPlot(sce3, "TNFRSF12A", split.by = "condition", cols = c(wt_col, ko_col),  pt.size = 0)
v3 <- VlnPlot(sce3, "JUN", split.by = "condition", cols = c(wt_col, ko_col),  pt.size = 0)
v1+v2+v3
dev.off()

#Figure 5SG 4E
Idents(sce1) <-sce1$condition 
markers1 <- FindAllMarkers(sce1, min.pct = 0.05, min.diff.pct = 0.02, logfc.threshold = .05, test.use = "MAST")
top <- markers1 %>% dplyr::filter(cluster=="WT")
keyvals <- ifelse(
  top$avg_log2FC < -.10, '#F2AD00',
  ifelse(top$avg_log2FC > .10, '#00A08A',
         'grey50'))
keyvals[is.na(keyvals)] <- 'grey50'
names(keyvals)[keyvals == '#00A08A'] <- 'Wildtype'
names(keyvals)[keyvals == 'grey50'] <- 'NS'
names(keyvals)[keyvals == '#F2AD00'] <- 'Knockout'
apotosis1 <-c("FAS",
              
              "TNFRSF10A",
              "TNFRSF10B",
              "TNFRSF10C",
              "TNFRSF10D",
              "TNFRSF11B",
              "TNFSF10",
              "TNFRSF1A",
              "TNFRSF12A",
              "FADD",
              "CFLAR",
              "CASP1",
              "CASP2",
              "CASP3",
              "CASP4",
              
              "CASP6",
              "CASP7",
              "CASP8",
              "CASP9",
              "CASP10",
              "CASP14",
              "NAIP",
              "BIRC2",
              "BIRC3",
              "XIAP",
              "BIRC5",
              "BIRC6",
              "BIRC7",
              "BCL2",
              "MCL1",
              "BCL2L1",
              "BCL2L2",
              "BCL2A1",
              "BCL2L10",
              "BAX",
              "BAK1",
              "BOK",
              "BID",
              "BCL2L11",
              "BMF",
              "BAD",
              "BIK",
              "HRK",
              "PMAIP1",
              "BNIP3",
              "BNIP3L",
              "BCL2L14",
              "BBC3",
              "BCL2L12",
              "BCL2L13",
              "APAF1",
              "CYCS",
              "DIABLO",
              "HTRA2",
              "AIFM1",
              "ENDOG",
              "CARD8",
              "CARD6",
              "NOX5",
              "TP53",
              "CDKN1A",
              "CDKN1B",
              "CDKN2A",
              "CDKN2B")
ofinterest <- c("MKI67",
                "TOP2A",
                "NES",
                "MAP2",
                "TH",
                "EN1",
                "NR4A2",
                "HIF1A",
                "TP53",
                "RELA",
                "NFKB1",
                "NFKB2",
                "JUN",
                "FOS", "LMO3", "HES5")
geneofinterest <- c(ofinterest, apotosis1)
#Label genes of interest, change Y axis, Save Size

samplabels1 <- c("PEG10", "DYNC1I2", "JUN",
                 
                 "NMT1",
                 "PPDPF",
                 "PHPT1",
                 "NDUFA13",
                 "NDUFC2",
                 "NDUFA11",
                 "NME2",
                 "NDUFA12",
                 "RPS27L",
                 "INSM1",
                 "FABP5",
                 "HSPA1A",
                 "WLS",
                 "MICOS13",
                 "ASCL1",
                 "SOX2")
pdf("Figure5sg_alt.pdf", width = 5, height = 5)
EnhancedVolcano(top,
                lab = top$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'All WT vs KO Top DE Genes',
                pCutoff = .0001,
                FCcutoff = .10,
                pointSize = 2.0,
                labSize = 3.3,
                col=c('blue', 'green', 'black', 'red3', "yellow", "brown", "grey"),
                colAlpha = .8,
                selectLab = samplabels1,
                colCustom = keyvals,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                lengthConnectors = unit(0.01, "npc"),
                cutoffLineType = "dotted", 
                cutoffLineCol = "#00000066") + xlim(-.4, .4)
                #titleLabSize = 10,
                #subtitleLabSize = 9,
                #captionLabSize = 9) 
               
dev.off()
sce4 <- subset(sce1, ident == 3 | ident == 5 | ident == 6)
Idents(sce4) <-sce4$condition

markers2 <- FindAllMarkers(sce4, min.pct = 0.01, min.diff.pct = 0.01, logfc.threshold = .05, test.use = "MAST")
top2 <- markers2 %>% dplyr::filter(cluster=="WT")
keyvals <- ifelse(
  top2$avg_log2FC > .10 & top2$p_val_adj < .01, '#00A08A',
  ifelse(top2$avg_log2FC < -.10 & top2$p_val_adj < .01, '#F2AD00',
         'grey50'))
keyvals[is.na(keyvals)] <- 'grey50'
names(keyvals)[keyvals == '#00A08A'] <- 'Wildtype'
names(keyvals)[keyvals == 'grey50'] <- 'NS'
names(keyvals)[keyvals == '#F2AD00'] <- 'Knockout'
samplabels<- c("RPS27L", "PHLDA3", "TUBB2A","TUBB2B","RPL27A","PHPT1","CDKN1A","BAX",
      "NMT1",
      "BBC3",
      "CDKN2B",
      "TNFRSF10B")
pdf("Figure4E_alt.pdf", width = 5, height = 5)
EnhancedVolcano(top2,
                lab = top2$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'WT vs KO Top DE Genes',
                pCutoff = .01,
                FCcutoff = .10,
                pointSize = 1.0,
                labSize = 3.3,
                selectLab = samplabels,
                col=c('blue', 'green', 'black', 'red3', "yellow", "brown", "grey"),
                colAlpha = .8,
                colCustom = keyvals,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                lengthConnectors = unit(0.01, "npc"),
                cutoffLineType = "dotted", 
                cutoffLineCol = "#00000066",
                #titleLabSize = 10,
                #subtitleLabSize = 8,
                #captionLabSize = 8
)  + xlim(-.75, .75) 
dev.off()
#Barplot and Gene plot 5SE 5SD
#FIXXXX TO BAR PLOT
table2 <- as.matrix(table(sce1@meta.data$condition, sce1@meta.data$cluster))
table3 <- data.table(t(t(table2) / colSums(table2)))
colnames(table3) <- c("Condition", "Cluster", "Fraction")
pdf("5SD_ALT.pdf", width = 3, height = 3)
#par(mar=c(4.1, 2.1, 4.1, 2.1))
#barplot(t(t(table2) / colSums(table2)), col=c(ko_col,wt_col), horiz=T, main="Proportion of Cluster by Group") +theme(text=element_text(size=16))
ggplot(data=table3, aes(x=Cluster, y=Fraction, fill=Condition)) +
  geom_bar(stat="identity") + scale_fill_manual(values = c(ko_col, wt_col))



dev.off()


sce1 <- sce2 # STOP HERE
Idents(sce1) <- sce1$condition
sce1 <- subset(x = sce1, idents = "KO") 

th <- length(WhichCells(sce1, slot = 'counts', expression = TH > 0)) /  length(WhichCells(sce1, slot = 'counts'))
mki67 <- length(WhichCells(sce1, slot = 'counts', expression = MKI67 > 0)) /  length(WhichCells(sce1, slot = 'counts'))
map2 <- length(WhichCells(sce1, slot = 'counts', expression = MAP2 > 0)) /  length(WhichCells(sce1, slot = 'counts'))
all_ko <- c(map2, th, mki67)
sce1 <- sce2
sce1 <- subset(x = sce1, idents = "WT")
th <- length(WhichCells(sce1, slot = 'counts', expression = TH > 0)) /  length(WhichCells(sce1, slot = 'counts'))
mki67 <- length(WhichCells(sce1, slot = 'counts', expression = MKI67 > 0)) /  length(WhichCells(sce1, slot = 'counts'))
map2 <- length(WhichCells(sce1, slot = 'counts', expression = MAP2 > 0)) /  length(WhichCells(sce1, slot = 'counts'))
all_wt <- c(map2, th, mki67)
temp <- data.table(
  Gene = c("MAP2", "TH", "MKI67"),
  Group =c("WT", "WT", "WT", "KO", "KO", "KO"),
  Values = c(all_wt, all_ko)
)

pdf("5SE.pdf", width = 4, height = 3)
ggplot(temp, aes(x = Gene, y = Values, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + ggtitle("Fraction of Cells Expressing Each Gene") + ggeasy::easy_center_title() +
  theme(plot.title = element_text(face="bold")) +theme(text=element_text(size=8)) + scale_fill_manual(values = c(ko_col, wt_col))
dev.off()

#FeaturePlots 5SF
sce1 <- sce2
sce1.list = list(sce1[, sce1@meta.data$condition == "WT"], sce1[, sce1@meta.data$condition == "KO"])

pdf("5SF.pdf",width = 7, height = 5)
#MAP2
f1 <-FeaturePlot(sce1.list[[1]], feature="MAP2", order=T, pt.size=1.0) + labs(title="MAP2 (WT)") + scale_colour_gradientn(colours = terrain.colors(10))
f2 <-FeaturePlot(sce1.list[[2]], feature="MAP2", order=T, pt.size=1.0) + labs(title="MAP2 (KO)") + scale_colour_gradientn(colours = terrain.colors(10))
#PBX1
f3 <-FeaturePlot(sce1.list[[1]], feature="PBX1", order=T, pt.size=1.0) + labs(title="PBX1 (WT)") + scale_colour_gradientn(colours =terrain.colors(10))
f4 <-FeaturePlot(sce1.list[[2]], feature="PBX1", order=T, pt.size=1.0) + labs(title="PBX1 (KO)") + scale_colour_gradientn(colours = terrain.colors(10))
#MKI67
f5 <-FeaturePlot(sce1.list[[1]], feature="MKI67", order=T, pt.size=1.0) + labs(title="MKI67 (WT)") + scale_colour_gradientn(colours = terrain.colors(10))
f6 <-FeaturePlot(sce1.list[[2]], feature="MKI67", order=T, pt.size=1.0) + labs(title="MKI67 (KO)") + scale_colour_gradientn(colours = terrain.colors(10))
(f1/f2) | (f3/f4) | (f5/f6)
dev.off()

#La Manno
combined_sce1 <- as.SingleCellExperiment(sce1)
combined_sce1$clusters <- sce1$clusters
combined_sce_lamanno <- scRNAseq::LaMannoBrainData("human-embryo")


# SingleR() expects reference datasets to be normalized and log-transformed.


combined_sce_lamanno <- logNormCounts(combined_sce_lamanno) 
predictions2 <- SingleR(test=combined_sce1, assay.type.test=2, 
                        ref=combined_sce_lamanno, labels =combined_sce_lamanno$Cell_type, de.method="wilcox", clusters=combined_sce1$clusters)
predictions2$pruned.labels #Didn't assign
sce5 <- as.Seurat(combined_sce1)
Idents(sce5) <- sce5$clusters
sce5 <- RenameIdents(object = sce5, `1` = "hNbM", `2` = "hNbM", `3` = "hProgFPL", '4'= "hNbM", "5" = "hProgFPL", "6" = "hProgFPL", "7" = "hPeric")
pdf("Figure4c.pdf", width = 4, height = 4)
DimPlot(sce5)
dev.off()

#4H by celltype
combined_sce_lamanno1  <- subset(combined_sce_lamanno , ,Cell_type == c('hNbM') | Cell_type == c('hNbML1') | Cell_type == c('hProgFPL'))
predictions3 <- SingleR(test=combined_sce1, assay.type.test=2, 
                        ref=combined_sce_lamanno1, labels =combined_sce_lamanno1$Cell_type, de.method="wilcox")
combined_sce1$pruned.labels <- predictions3$pruned.labels
sce6 <- as.Seurat(combined_sce1)

pdf("Figure4H.pdf", width = 4, height = 6.5)
#cols2 <- wes_palette("Rushmore1",3)
v1 <- VlnPlot(sce6, "BAX", group.by = "pruned.labels", cols = c(wt_col, ko_col), pt.size = 0, split.by = "condition")
v2 <- VlnPlot(sce6, "TNFRSF12A", group.by = "pruned.labels",  cols = c(wt_col, ko_col),  pt.size = 0, split.by = "condition")
v3 <- VlnPlot(sce6, "JUN", group.by = "pruned.labels", cols = c(wt_col, ko_col),  pt.size = 0, split.by = "condition")
v1 + v2 + v3
dev.off()
#CellLine
#Load Cellline
cellline <- readRDS("cellline_clean.rds")


pdf("Figure4G_alt.pdf",width = 6, height = 4)
#HES5
f1 <- FeaturePlot(sce1, feature="HES5", order=T, pt.size=1.0) + labs(title="HES5") + scale_colour_gradientn(colours = terrain.colors(10))
f2 <- FeaturePlot(cellline, feature="HES5", order=T, pt.size=1.0) + labs(title="HES5") + scale_colour_gradientn(colours = terrain.colors(10))
f1|f2
dev.off()

#ALTERNATE FIGURE REQUESTS
sce4 <- subset(sce1, ident == 1 |ident == 2 | ident == 4)
Idents(sce4) <-sce4$condition

markers2 <- FindAllMarkers(sce4, min.pct = 0.01, min.diff.pct = 0.01, logfc.threshold = .05, test.use = "MAST")
top2 <- markers2 %>% dplyr::filter(cluster=="WT")
keyvals <- ifelse(
  top2$avg_log2FC >= .10 & top2$p_val_adj < .01, '#00A08A',
  ifelse(top2$avg_log2FC <= -.10 & top2$p_val_adj < .01, '#F2AD00',
         'grey50'))
keyvals[is.na(keyvals)] <- 'grey50'
names(keyvals)[keyvals == '#00A08A'] <- 'Wildtype'
names(keyvals)[keyvals == 'grey50'] <- 'NS'
names(keyvals)[keyvals == '#F2AD00'] <- 'Knockout'
samplabels<- c("PEG10",
               "NMT1",
               "HSPA8", "RPS27L",
               "NDUFA12",
               "DYNLT1",
               "PDCD5",
               "PHPT1",
               "KRTCAP2", "SYT14",
               "SOX2",
               "HSPA1B",
               "HSPA1A","HSPA2A","DNAJB1","HSPA8")
pdf("Figures_volcano2,4lg2fc.10_.pdf", width = 5, height = 5)
EnhancedVolcano(top2,
                lab = top2$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'WT vs KO Top DE Genes',
                pCutoff = .01,
                FCcutoff = .10,
                pointSize = 1.0,
                labSize = 3.3,
                selectLab = samplabels,
                col=c('blue', 'green', 'black', 'red3', "yellow", "brown", "grey"),
                colAlpha = .8,
                colCustom = keyvals,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                lengthConnectors = unit(0.01, "npc"),
                cutoffLineType = "dotted", 
                cutoffLineCol = "#00000066",
                #titleLabSize = 10,
                #subtitleLabSize = 8,
                #captionLabSize = 8
)  + xlim(-.8, .8) 
dev.off()


##### ALTERNATE VIOLIN PLOTS
sce3 <- subset(sce1, ident = 7, invert = TRUE)
pdf("FigureALTVIOLIN8-9.pdf", width = 3, height = 6) #ATF3, ATF4, JUN, FOS, HSPD1, HSPA1A, HSPA2A, and DNAJB1  HSPA8 	DNAJB2 ATF2
v1 <- VlnPlot(sce3, "ATF3", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col), pt.size = 0)
v2 <- VlnPlot(sce3, "ATF4", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v3 <- VlnPlot(sce3, "JUN", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v4 <- VlnPlot(sce3, "FOS", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col), pt.size = 0)
v5 <- VlnPlot(sce3, "HSPD1", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v6 <- VlnPlot(sce3, "HSPA1A", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v8 <- VlnPlot(sce3, "DNAJB1", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v9 <- VlnPlot(sce3, "HSPA8", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v10 <- VlnPlot(sce3, "DNAJB2", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v11 <- VlnPlot(sce3, "ATF2", split.by = "condition", split.plot = TRUE, cols = c(wt_col, ko_col),  pt.size = 0)
v8+v9
dev.off()
