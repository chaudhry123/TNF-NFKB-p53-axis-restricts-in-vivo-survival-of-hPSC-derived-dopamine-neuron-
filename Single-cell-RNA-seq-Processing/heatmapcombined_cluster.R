#Import libraries
rm(list = ls())
library(DropletUtils)
library(scater)
#library(AnnotationHub)
library(scran)
#library(DropletTestFiles)
#library(EnsDb.Hsapiens.v86)
#library(pheatmap)
library(ggplot2)
library(Seurat)
#library(intrinsicDimension)
#library(batchelor)
library(dplyr)
#library(animation)
library(viridis)
library(EnhancedVolcano)
library(formattable)
library(pheatmap)
library(ape)
library(ComplexHeatmap)
#("Scillus")
addTaskCallback(functilibraryon(...) {set.seed(100);TRUE}) #Permanently sets seed so is reproducible

#ATTENTION YOU SUBBED IN GENE SYMBOLS FOR ROWNAMES AFTER MAKING UNIQUE
#Markers
mDA.mature.markers <- c("TH", "NR4A2", "PITX3", "CHRNA4", "SNCA", "PBX1", "SLC6A3", "SLC18A2", "SNCG", "ADCYAP1")
mDA.precursor.markers <- c("EN1", "FOXA1","FOXA2","LMX1A","LMX1B","SOX6","WNT1", "WNT5A",
                           "POSTN", "RELN", "SPOCK1")
A9.markers <- c("ALDH1A1", "LMO3", "SLC17A6", "NDNF", "SATB1","KCNJ6")
A10.markers <- c("OTX2", "CALB1", "CALB2", "CCK", "VGF", "VIP")
non.mDA.markers <- c("NKX2-1", "NKX2-2","HOXB2","GBX2", "PAX6", "FOXG1",
                     "DBX1", "BARHL1","BARHL2","PITX2","DBH", "POU5F1",
                     "NANOG","DLX2","ISL1","GATA3","SLC6A4","SIM1",
                     "POU4F1","TTR","COL1A1","SIX1","ACTA2","GFAP",
                     "ACTB", "EN2","TP53", "BBC3","CDKN1A", "ATOH1")
apotosis <-c("FAS",
             "FASLG",
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
             "CASP5",
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
                
                "RELA",
                "NFKB1",
                "NFKB2",
                "JUN",
                "FOS")
genesofinterest <- c(ofinterest, mDA.mature.markers, mDA.precursor.markers, apotosis, A9.markers, A10.markers, "PHPT1")
#Load Object
#combined_sce1 <- readRDS("umapfinalfinal.rds")

#Basic Clustering
#g <- buildSNNGraph(sce1, k=40, use.dimred = 'PCA')
#clust <- igraph::cluster_walktrap(g)$membership
#colLabels(sce1) <- factor(clust)

#table(colLabels(combined_sce1))
#plotUMAP(sce1, colour_by="label")
#Analysis in Seurat
#sampnames <- combined_sce2$Sample
#sampnames2 <- substr(sampnames, 0, 2)
combined_sce2 <- combined_sce1 

gene <- rowData(combined_sce2)$Symbol
temp <- make.unique(gene, sep ="_")
rownames(combined_sce2) <- temp
sce2 <- as.Seurat(combined_sce2)
tmpIdent <- combined_sce2$cluster
Idents(sce2) <- tmpIdent
sce1 <- sce2


#########
#VlnPlot(sce1, features = apotofinterest[14:26])

#Apotosis genes
#VlnPlot(sce1, features = )
#Dopamine Neurons
#VlnPlot(sce1, features = c("MKI67", "TOP2A")) #Look into MKI67 and TOP2A
#VlnPlot(sce1, features = ofinterest[13:14])
#Of interest
#FeaturePlot(sce1, features = ofinterest[7:15])
#Of interest

##############
#MAKE HEAT MAP, CLUSTERPLOTS, VLN
Idents(fayzan.4) <- fayzan.4$cluster
markers <- FindAllMarkers(sce1, min.pct = 0.05, min.diff.pct = 0.05, logfc.threshold=.10, only.pos = TRUE)
top20 <- markers %>%  group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
sce1 <- ScaleData(object = fayzan.4, block.size=100)
#png("heatmaptopgroup.png", width = 800, height = 800)
mat<- temp3[unique(top20$gene), ] %>% as.matrix()
#mat<- t(scale(t(mat)))
cluster_anno<- Idents(sce1)
col_fun = circlize::colorRamp2(c(-1, 0, 2), c("#140b34" , "#84206b", "#f6d746"))
Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        row_dend_width = unit(20, "mm"),
        col = col_fun,
        row_names_gp = gpar(fontsize = 4),
        column_title_rot = 90,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)
dev.off()
#png("precursorheatmap.png", width = 800, height = 800)
#DoHeatmap(object = sce1, features = mDA.precursor.markers, label = TRUE, slot="scale.data")
#dev.off()
png("precursorviolin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = "mDA.precursor.markers"BAX, pt.size = 0, group.by = "cells")
dev.off()
png("matureviolin2group.png", width = 800, height = 800)
VlnPlot(sce2, features = "TH" , pt.size = 0)
dev.off()
png("a9violin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = A9.markers, pt.size = 0)
dev.off()
png("a10violin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = A10.markers, pt.size = 0)
dev.off()
png("ofinterestviolin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = sce1@meta.data$Sample, pt.size = 0)
dev.off()
#VlnPlot(sce1, features = apotosis[60:70])
#FeaturePlot(sce1, features = mDA.precursor.markers)
#FeaturePlot(sce1, features = mDA.mature.markers)
#FeaturePlot(sce1, features = A9.markers)
#FeaturePlot(sce1, features = A10.markers)
png("ofinterestpart1group.png", width = 800, height = 800)
FeaturePlot(sce1, features = "BAX")
dev.off()
#png("ofinterestpart2group.png", width = 800, height = 800)
FeaturePlot(sce1, features = apotosis)
dev.off()
#cluster.markers <- FindAllMarkers(sce1, features = apotosis)
png("apotosisheatmapgroup.png", width = 800, height = 800)
DoHeatmap(object = sce1, features = apotosis, label = TRUE, slot="scale.data") + scale_fill_gradientn(colors=inferno(256))
dev.off()
png("ofinterestheatmapgroup.png", width = 800, height = 800)
DoHeatmap(object = sce1, features = ofinterest, label = TRUE, slot="scale.data") + scale_fill_gradientn(colors=inferno(256))
dev.off()


#EXTRA
#unique(top11$gene)
all.markers <- unique(c(ofinterest, apotosis))#, #ofinterest)) ,apotosis))
markers <- FindAllMarkers(sce1,  min.pct = 0.01, min.diff.pct = 0.01, logfc.threshold = .01, only.pos = TRUE, features = genes5)
top14 <- markers %>% arrange(p_val_adj) #check that adjusted pvals are correct
top11 <- top14[1:116,]
#sce1 <- ScaleData(object = sce1, block.size=100)
#png("heatmaptopapopsamp.png", width = 800, height = 800)
#DoHeatmap(object = sce1, features = top10$gene, label = TRUE, slot="scale.data") + scale_fill_gradientn(colors=inferno(256))
#dev.off()
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"
customRed0 = "#994C4C"
DotPlot(sce1, features=c(unique(genes5)),  dot.scale = 6, cols =c("RdGy"))  + RotatedAxis() +theme(axis.text.x = element_text(color = "black", size = 9))
formattable(top11, align =c("l","c","c","c","c", "c", "r"), list(
  `Indicator Name` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `p_val_adj` = color_tile(customRed0,customRed), 
  `avg_log2FC` = color_tile(customGreen0, customGreen )
))

VlnPlot(sce1, features = c(genes5), pt.size = 0)
FeaturePlot(sce1, features =c("ALDH1A1", "LMO3"))
genes5 <- c("OTX2", "CALB1","LMO3", "ALDH1A1", "SLC182A", "SLC182A", "TH", "NR4A2", "PBX1" ,"ADCYAP1", "MAP2", "FOXA2", "JUN", "TP53",  "TNFRSF12A", "PUMA", "TNFRSF1A", "TNFRSF10D", "BAX", "TNFRSF10B", "TNFRSF11B", "MKI67", "HES5")

EnhancedVolcano(top10,
                lab = top10$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Top DE Genes',
                pCutoff = .01,
                FCcutoff = .10,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('blue', 'green', 'black', 'red3'),
                colAlpha = 1)


markers <- FindAllMarkers(sce3, min.diff.pct = 0.01, logfc.threshold=.01, min.pct = .02)
top10 <- markers %>% dplyr::filter(cluster=="hNbM")
top10 <- markers %>% dplyr::filter(p_val_adj < .05) %>%top_n(50, avg_log2FC)
sce1 <- ScaleData(object = sce1, block.size=100)
#png("heatmaptopgroup.png", width = 800, height = 800)
pheatmap(object = sce1, features = top10$gene, slim.col.label = TRUE, remove.key = TRUE) + scale_fill_gradientn(colors=inferno(256))
#dev.off()
BuildClusterTree(sce1)
PlotClusterTree(sce1)

#
library(org.Hs.eg.db)
markers <- markers[,7:6]
dfsample <- split(markers$gene,markers$cluster)
length(dfsample)

dfsample$`WT` = bitr(dfsample$`WT`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`KO` = bitr(dfsample$`KO`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#dfsample$`6` = bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



genelist <- list("WT" = dfsample$`WT`$ENTREZID, 
                 "KO" = dfsample$`KO`$ENTREZID
                # "3" = dfsample$`2`$ENTREZID,
                 #"4" = dfsample$`2`$ENTREZID,
                # "5" = dfsample$`2`$ENTREZID,
                # "6" = dfsample$`2`$ENTREZID
                 )
clusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(clusterplot)       
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = org.Hs.eg.db)
dotplot(GOclusterplot)    




pbx1 <- length(WhichCells(sce1, slot = 'counts', expression = MAP2 > 0 & PBX1 > 0)) /  length(WhichCells(sce1, slot = 'counts'))




#Segment on cluster 6 then on bax
#th <- WhichCells(sce1, slot = 'counts', expression = BAX > 0)
sce3 <- subset(sce1,  pruned_labels == "hNbM")
#Idents(sce6) <- tmpIdent
table(Idents(sce6))
#24859 18577 




#
library(scales) 
library(dplyr)
library(ggplot2)
library(tidyr)
#
TH_WT <- subset(fayzan.4, TH > 0 & fayzan.4$condition =="WT")
TH_KO <- subset(fayzan.4, TH > 0 & fayzan.4$condition =="KO")
MAP2_WT <- subset(fayzan.4, MAP2 > 0 & fayzan.4$condition =="WT")
MAP2_KO <- subset(fayzan.4, MAP2 > 0 & fayzan.4$condition =="KO")
MKI67_WT <- subset(fayzan.4, MKI67 > 0 & fayzan.4$condition =="WT")
MKI67_KO <-subset(fayzan.4, MKI67 > 0 & fayzan.4$condition =="KO")

TH_WT <- dim(TH_WT)[2]
TH_KO <- dim(TH_KO)[2]
MAP2_WT <- dim(MAP2_WT)[2]
MAP2_KO <- dim(MAP2_KO)[2]
MKI67_WT <- dim(MKI67_WT)[2]
MKI67_KO <-dim(MKI67_KO)[2]
#
fayzan.4$cluster <- Idents(fayzan.4)
meta <- fayzan.4@meta.data
d2 <- meta %>% 
  group_by(cluster, as.factor(condition)) %>%  #Cluster, sample
  summarise(count = n())

a <- table(meta$condition)
b <- c(a,a,a,a,a,a,a)
c <- d2$count
d <- c/b 
d2$perc <- d
d2$condition <- d2$`as.factor(condition)`


png("gene_Percent.png")
ggplot(a, aes(x = factor(gene), y = perc)) + 
  geom_bar(aes(y = perc, fill = factor(condition)), stat="identity", position = "dodge") + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Percentages") + guides(fill=guide_legend(title="Condition")) + xlab("Gene") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
dev.off()


all <- rbind(all_wt1,all_wt2, all_ko1, all_ko2)
rownames(all) <- c("WT1", "WT2", "KO1", "KO2")
data_base <- reshape(all,                        # Modify data for Base R barplot
                     idvar = "subgroup",
                     timevar = "group",
                     direction = "wide")

  
  
  
  
#New UMAP
VlnPlot(sce1, features = A10.markers, pt.size = 0, group.by = "pruned_labels", split.by = 'condition') 
VlnPlot(sce1, features = "JUN", pt.size = 0, group.by = "cells")
VlnPlot(sce1, features = "TNFRSF12A", pt.size = 0.1, group.by = "cells")
sce1$cluster <- Idents(sce1)
Idents(sce1) <- sce1$condition
sce4 <- subset(sce1, cluster == '7')

markers_7 <- FindMarkers(sce1, logfc.threshold= .05, min.pct = 0.05, min.diff.pct = 0.02, ident.1 ="WT", ident.2="KO", test.use = "MAST")

markers_1plus <- markers_1 %>% dplyr::filter(avg_log2FC > .05) %>% dplyr::filter(p_val < .01)

EnhancedVolcano(markers_1,
                lab = rownames(markers_1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #selectLab = genesofinterest,
                title = 'Top DE Genes for All Clusters',
                pCutoff = .01,
                FCcutoff = .10,
                pointSize = 3.0,
                labSize = 4.0,
                col=c('blue', 'green', 'black', 'red3'),
                colAlpha = 1,
                xlim =c(-.75,.75))
  
#UPSETR
sce4 <- subset(fayzan.4)
Idents(sce4) <- sce4$condition
markers_1 <- FindMarkers(sce4, logfc.threshold= .05, min.pct = 0.05, min.diff.pct = 0.02, ident.1 ="WT", ident.2="KO", test.use = "MAST")
markers_3_vol <- markers_3 %>% dplyr::filter(avg_log2FC > .05) %>% dplyr::filter(p_val < .01)
markers_1neg <- markers_1 %>% dplyr::filter(avg_log2FC < -.05) %>% dplyr::filter(p_val < .01)
a <- list(one = c(rownames(markers_1plus)), two = c(rownames(markers_2plus)), three = c(rownames(markers_3plus)), four = c(rownames(markers_4plus)),five = c(rownames(markers_5plus)),six = c(rownames(markers_6plus)),seven = c(rownames(markers_7plus)))
b <- list(one = c(rownames(markers_1neg)), two = c(rownames(markers_2neg)), three = c(rownames(markers_3neg)), four = c(rownames(markers_4neg)),five = c(rownames(markers_5neg)),six = c(rownames(markers_6neg)),seven = c(rownames(markers_7neg)))
upset(fromList(a), order.by = "freq", nsets=7)
upset(fromList(b), order.by = "freq", nsets=7)


#Plot
sce1 <- readRDS("cellline_clean.rds")
fayzan.5 <- subset(fayzan.4, cluster != '7')
VlnPlot(sce4, features = c("BAX", "TNFRSF12A"), pt.size = 0, split.by = "condition")
fayzan.4.list = list(fayzan.4[, fayzan.4@meta.data$condition == "WT"], fayzan.4[, fayzan.4@meta.data$condition == "KO"])

pdf("FeaturePlot.infernocolor.20220512_NPC2.pdf", width=7.5, height=7.5)
# TP53
FeaturePlot(sce1, feature="TNFRSF12A", order=T, pt.size=3.0) + labs(title="TNFRSF12A") + scale_color_gradientn(colors=c("grey95", rev(magma(99))), limits=c(0, 5)) 

FeaturePlot(fayzan.4.list[[1]], feature="NPC2", order=T, pt.size=3.0) + labs(title="NPC2 (WT)") + scale_color_gradientn(colors=c("grey95", rev(magma(99))), limits=c(0, 5)) 

FeaturePlot(fayzan.4.list[[2]], feature="NPC2", order=T, pt.size=3.0) + labs(title="NPC2 (KO)") + scale_color_gradientn(colors=c("grey95", rev(magma(99))), limits=c(0, 5))
dev.off()

#Updated Violin
pdf("Panel-D.VlnPlot.APOPGENES_ALL.20220425.pdf", width=7, height=7)
p1 = VlnPlot(fayzan.4, feature=c("BAX", "TNFRSF12A", "BBC3", "PHPT1", "TNFRSF10A"), split.by="condition", split.plot=T, cols=c("green", "pink"), pt.size=0)
p1 + geom_jitter(mapping=aes(color=split), data=p1$data, size=0.05, position=position_jitterdodge(jitter.width=0.4, dodge.width=0.9)) + guides(color="none") + scale_color_manual(breaks=c("WT", "KO"), values=c("green", "pink"))
dev.off()


tnf <- 
