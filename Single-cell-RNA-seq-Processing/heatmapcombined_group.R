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
#("Scillus")
addTaskCallback(functilibraryon(...) {set.seed(100);TRUE}) #Permanently sets seed so is reproducible

#ATTENTION YOU SUBBED IN GENE SYMBOLS FOR ROWNAMES AFTER MAKING UNIQUE
#Markers
mDA.mature.markers <- c("TH", "NR4A2", "PITX3",  "SNCA", "PBX1",  "SLC18A2",  "ADCYAP1")
mDA.precursor.markers <- c("EN1", "FOXA1","FOXA2","LMX1A","LMX1B","SOX6","WNT1", "WNT5A",
                           "POSTN", "SPOCK1") #"CHRNA4", "SLC6A3", "SNCG", "RELN",
A9.markers <- c("LMO3",  "NDNF", "SATB1","KCNJ6") #"ALDH1A1") "SLC17A6",
A10.markers <- c("OTX2", "CALB1", "CALB2", "CCK", "VGF")# "VIP")
non.mDA.markers <- c("NKX2-1",  "PAX6", "FOXG1",
                      "BARHL2","PITX2", "POU5F1",
                     "DLX2","ISL1","SLC6A4",
                     "POU4F1","TTR","COL1A1","SIX1","ACTA2",
                     "ACTB", "TP53", "BBC3","CDKN1A") #"EN2","GFAP","GATA3", "NANOG", "BARHL1","DBH","SIM1", ,"HOXB2" "NKX2-2", "GBX2", "DBX1",
apotosis <-c("FAS",
             #"FASLG",
             "TNFRSF10A",
             "TNFRSF10B",
           #  "TNFRSF10C",
             "TNFRSF10D",
             "TNFRSF11B",
             #"TNFSF10",
             "TNFRSF1A",
             "TNFRSF12A",
             "FADD",
             "CFLAR",
        #"CASP1",
             "CASP2",
             "CASP3",
           # "CASP4",
          #   "CASP5",
             "CASP6",
             "CASP7",
            # "CASP8",
             "CASP9",
            # "CASP10",
            # "CASP14",
             "NAIP",
             "BIRC2",
             "BIRC3",
             "XIAP",
             "BIRC5",
             "BIRC6",
            # "BIRC7",
             "BCL2",
             "MCL1",
             "BCL2L1",
             "BCL2L2",
          #   "BCL2A1",
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
            # "CDKN2A",
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
#Load Object
#sce1 <- readRDS("processed_file")

#Basic Clustering
#g <- buildSNNGraph(sce1, k=40, use.dimred = 'PCA')
#clust <- igraph::cluster_walktrap(g)$membership
#colLabels(sce1) <- factor(clust)

#table(colLabels(combined_sce1))
#plotUMAP(sce1, colour_by="label")
#Analysis in Seurat
sampnames <- combined_sce1$Sample
sampnames2 <- substr(sampnames, 0, 2)
combined_sce2 <- combined_sce1

gene <- rowData(combined_sce2)$Symbol
temp <- make.unique(gene, sep ="_")
rownames(combined_sce2) <- temp
sce2 <- as.Seurat(combined_sce2)
tmpIdent <- sampnames2
Idents(sce2) <- tmpIdent
#
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
markers <- FindAllMarkers(sce1, only.pos=TRUE)
top10 <- markers %>% dplyr::group_by(cluster) %>% top_n(10)
sce1 <- ScaleData(object = sce1, block.size=100)
png("heatmaptopgroup.png", width = 800, height = 800)
DoHeatmap(object = sce1, features = markers$gene, label = TRUE, slot="scale.data") + scale_fill_gradientn(colors=inferno(256))
dev.off()
#png("precursorheatmap.png", width = 800, height = 800)
#DoHeatmap(object = sce1, features = mDA.precursor.markers, label = TRUE, slot="scale.data")
#dev.off()
png("precursorviolin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = mDA.precursor.markers, pt.size = 0)
dev.off()
png("matureviolin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = mDA.mature.markers, pt.size = 0)
dev.off()
png("a9violin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = A9.markers, pt.size = 0)
dev.off()
png("a10violin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = A10.markers, pt.size = 0)
dev.off()
png("ofinterestviolin2group.png", width = 800, height = 800)
VlnPlot(sce1, features = all.markers, pt.size = 0)
dev.off()
#VlnPlot(sce1, features = apotosis[60:70])
#FeaturePlot(sce1, features = mDA.precursor.markers)
#FeaturePlot(sce1, features = mDA.mature.markers)
#FeaturePlot(sce1, features = A9.markers)
#FeaturePlot(sce1, features = A10.markers)
png("ofinterestpart1group.png", width = 800, height = 800)
FeaturePlot(sce1, features = ofinterest[0:7])
dev.off()
#png("ofinterestpart2group.png", width = 800, height = 800)
FeaturePlot(sce1, features = all.markers)
sce1dev.off()
#cluster.markers <- FindAllMarkers(sce1, features = apotosis)
png("apotosisheatmapgroup.png", width = 800, height = 800)
DoHeatmap(object = sce1, features = apotosis, label = TRUE, slot="counts") + scale_fill_gradientn(colors=inferno(256))
dev.off()
png("ofinterestheatmapgroup.png", width = 800, height = 800)
DoHeatmap(object = sce1, features = ofinterest, label = TRUE, slot="counts") + scale_fill_gradientn(colors=inferno(256))
dev.off()


#EXTRA
#unique(top11$gene)
all.markers <- unique(c(ofinterest, apotosis, mDA.mature.markers, mDA.precursor.markers, A10.markers, A9.markers, non.mDA.markers))#, #ofinterest)) ,apotosis))
markers <- FindAllMarkers(sce1, features = unique(all.markers), min.pct = 0.01, min.diff.pct = 0.01, logfc.threshold = .01)
top14 <- markers %>% dplyr::filter(cluster=="WT1")   #check that adjusted pvals are correct
#top11 <- top14[1:116,]
#sce1 <- ScaleData(object = sce1, block.size=100)
#png("heatmaptopapopsamp.png", width = 800, height = 800)
#DoHeatmap(object = sce1, features = top10$gene, label = TRUE, slot="scale.data") + scale_fill_gradientn(colors=inferno(256))
#dev.off()
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"
customRed0 = "#994C4C"
DotPlot(sce1, features=c(unique(top14$gene)),  dot.scale = 6, cols =c("RdGy"))  + RotatedAxis() +theme(axis.text.x = element_text(color = "black", size = 7))
formattable(top14, align =c("l","c","c","c","c", "c", "r"), list(
  `Indicator Name` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `p_val_adj` = color_tile(customRed0,customRed), 
  `avg_log2FC` = color_tile(customGreen0, customGreen )
))

VlnPlot(sce1, features = c("TNFRSF12A"), pt.size = 0)
FeaturePlot(sce1, features =c("TNFRSF12A"))
keyvals.color <- top14$cluster
names(keyvals.color)[keyvals.color == 1] <- 'Cluster 1'
names(keyvals.color)[keyvals.color == 2] <- 'Cluster 2'
names(keyvals.color)[keyvals.color == 3] <- 'Cluster 3'
names(keyvals.color)[keyvals.color == 4] <- 'Cluster 4'
names(keyvals.color)[keyvals.color == 5] <- 'Cluster 5'
names(keyvals.color)[keyvals.color == 6] <- 'Cluster 6'
names(keyvals.color)[keyvals.color == 7] <- 'Cluster 7'
#top14$p_val_adj <- top14$p_val
EnhancedVolcano(top14,
                lab = top14$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'WT1 vs KO1 Top DE Genes from Genes of Interest',
                pCutoff = .01,
                FCcutoff = .05,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('blue', 'green', 'black', 'red3', "yellow", "brown", "grey"),
                colAlpha = 1)
                #colCustom = keyvals.color)

markers1 <- FindAllMarkers(sce1, features= c(ofinterest, apotosis), min.pct = 0.01, min.diff.pct = 0.01, logfc.threshold = .02)
top <- top14 %>% dplyr::select(ends_with())
newdata <- markers[cluster==1, ]
top <- markers %>% dplyr::filter(cluster==1)
sce1 <- ScaleData(object = sce1, block.size=100)
#png("heatmaptopgroup.png", width = 800, height = 800)
pheatmap(object = sce1, features = top10$gene, slim.col.label = TRUE, remove.key = TRUE) + scale_fill_gradientn(colors=inferno(256))
#dev.off()
BuildClusterTree(sce1)
PlotClusterTree(sce1)


#

markers <- markers1[,7:6]
dfsample <- split(markers$gene,markers$cluster)
length(dfsample)

dfsample$`Graft` = bitr(dfsample$`Graft`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Cellline` = bitr(dfsample$`Cellline`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


genelist <- list("Graft" = dfsample$`Graft`$ENTREZID, 
                 "Cellline" = dfsample$`Cellline`$ENTREZID)
clusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(clusterplot)       
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = org.Hs.eg.db)
dotplot(GOclusterplot)    


tnfgenes <- c("A2M",
  "ABCB1",
  "ABHD3",
  "ABTB1",
  "AC090186.1",
  "ADGRV1",
  "ADORA2B",
  "ADPRM",
  "ADRB2",
  "AEBP1",
  "AIFM3",
  "ANG",
  "ANGPT2",
  "ANKRD37",
  "AP1G2",
  "AP1S2",
  "APOA1",
  "AQP4",
  "ARHGAP29",
  "ARHGAP36",
  "ARL13B",
  "ARL4A",
  "ARL4D",
  "ARL6",
  "ARSK",
  "ATF1",
  "ATG14",
  "ATG16L2",
  "ATXN3",
  "B2M",
  "BAG3",
  "BCL3",
  "BGN",
  "BIK",
  "BMP5",
  "BROX",
  "BTG1",
  "BTG3",
  "C7orf33",
  "C7orf61",
  "CAPN9",
  "CAV1",
  "CCDC102B",
  "CCDC186",
  "CCDC28A",
  "CCDC66",
  "CCL2",
  "CCL20",
  "CD44",
  "CD70",
  "CD9",
  "CDH9",
  "CDK2AP2",
  "CDKN2A",
  "CENPC",
  "CEP162",
  "CFC1",
  "CHAC1",
  "CHORDC1",
  "CIART",
  "CISD2",
  "CLEC1A",
  "CLK1",
  "COL21A1",
  "CRH",
  "CRYAB",
  "CSMD3",
  "CTHRC1",
  "CTSC",
  "CXCL2",
  "CXCL5",
  "CXCL6",
  "CXCR4",
  "DDIT3",
  "DEPDC7",
  "DEPP1",
  "DHH",
  "DKK1",
  "DLC1",
  "DNAJB1",
  "DNAJC3",
  "DRD4",
  "DSEL",
  "ELF1",
  "ELL3",
  "EMP1",
  "EMP2",
  "EMP3",
  "ENO3",
  "EPCAM",
  "EPS8L2",
  "ERO1B",
  "ESCO1",
  "F2RL1",
  "FBXO32",
  "FHL2",
  "FOXG1",
  "FST",
  "FSTL5",
  "FXYD5",
  "GAD2",
  "GAL",
  #"GBX2",
  "GPR37",
  "GPR4",
  "GRIA2",
  "GRIA3",
  "GSDMB",
  "H1-2",
  "H2AC6",
  "H2BC5",
  "HBA2",
  "HERPUD1",
  "HES5",
  "HEXIM1",
  "HEY1",
  "HEY2",
  "HEYL",
  "HIVEP1",
  "HIVEP3",
  "HSD17B11",
  "HSD17B14",
  "HSD17B2",
  "HSPA1A",
  "HSPA1B",
  "HSPA6",
  "HSPB1",
  "HSPG2",
  "HSPH1",
  "ID4",
  "IFITM3",
  "IGF1",
  "IGFBP6",
  "IL13RA2",
  "IL32",
  "ING3",
  "IRS2",
  "KCTD4",
  "KCTD6",
  "KLF3",
  "KLF4",
  "KLF5",
  "KRAS",
  "LAMTOR3",
  "LAYN",
  "LDHA",
  "LGALS1",
  "LGALS3",
  "LGI1",
  "LRRIQ3",
  "LURAP1L",
  "LYG1",
  "MALRD1",
  "MAMDC4",
  "MANF",
  "MC5R",
  "MFGE8",
  "MMP7",
  "MMP9",
  "MRPL42",
  "MT-CO1",
  "MTRNR2L8",
  "NECTIN3",
  "NEIL1",
  "NPB",
  "NPC2",
  "NR4A1",
  "NTAN1",
  "NTF3",
  "NTS",
  "OSER1",
  "OXTR",
  "PARPBP",
  "PCDHGA12",
  "PCMTD1",
  "PGAP1",
  "PGF",
  "PHLDA2",
  "PKD1L3",
  "PLA2G4A",
  "PLAAT3",
  "PLAUR",
  "PLIN2",
  "PLOD2",
  "PLS1",
  "PMCH",
  "PMEL",
  "PNRC1",
  "POSTN",
  "PPP1R13L",
  "PRLR",
  "PTCHD1",
  "PTPN12",
  "QRFPR",
  "RAB25",
  "RABGGTB",
  "RBP1",
  "RBPMS",
  "RELB",
  "RENBP",
  "RGS2",
  "RIT1",
  "RIT2",
  "RNF113B",
  "RNF213",
  "RP13-314C10.6",
  "RPL22L1",
  "RPS10",
  "RRAD",
  "RSRP1",
  "RTKN2",
  "S100A10",
  "S100A11",
  "S100A13",
  "S100A16",
  "S100A2",
  "S100A3",
  "SAAL1",
  "SAT1",
  "SCGN",
  "SCRG1",
  "SDF2L1",
  "SDHD",
  "SERPINA3",
  "SERTAD1",
  "SESN3",
  "SETD4",
  "SFR1",
  "SHISAL2B",
  "SIGIRR",
  "SLC16A3",
  "SLC35G5",
  "SLC44A3",
  "SLC49A3",
  "SNAI1",
  "SP110",
  "SPACA6",
  "SPAG4",
  "SPCS3",
  "SUMO4",
  "TAF1D",
  "TCEAL9",
  "TCP11L2",
  "TDRD7",
  "TENT5A",
  "TEP1",
  "TFF1",
  "TFPI",
  "TGFB1",
  "TIFA",
  "TM4SF1",
  "TM4SF18",
  "TMEM156",
  "TMEM263",
  "TMEM33",
  "TMEM74",
  "TMTC2",
  "TNFRSF4",
  "TNFRSF6B",
  "TNFRSF9",
  "TNFSF12",
  "TNFSF9",
  "TRMO",
  "TSPO",
  "TUBA1C",
  "TYW1B",
  "UBD",
  "USPL1",
  "VAMP4",
  "VCAM1",
  "XBP1",
  "YPEL3",
  "ZAN",
  "ZC3H12A",
  "ZCCHC10",
  "ZFAND2A",
  "ZMYM5",
  "ZNF235",
  "ZNF490")

tnfgenes1 <- tnfgenes[tnfgenes %in% gene]

#GO ANALYSIS PART 2
markers1 <- FindAllMarkers(sce1, min.pct = 0.01, min.diff.pct = 0.01, logfc.threshold = .01)

top <- markers1 %>% dplyr::filter(cluster=="KO") #%>% filter(p_val_adj < .05)
markers <- top1[,c(7,2, 5)]
               
RA_output <- run_pathfindR(markers, gene_sets = "GO-All")
rtrial <- RA_output %>% 
  arrange(desc(Fold_Enrichment)) %>% slice(1:40)
RA_selected <- subset(RA_output, ID %in% go)
enrichment_chart(RA_selected)
#saveRDS(RA_output, "wtvscelllinego.rds")
go <- c("GO:0042981", "GO:0097193", "GO:2001244", "GO:0042771", "GO:0031072", "GO:0000045")
go <- c("GO:0002039", "GO:0031072", "GO:0005776", "GO:0007249", "GO:0051059",  "GO:0042981", "GO:0006915", "GO:0043065", "GO:0007250", "GO:2001238")
go <- c("GO:0002039", "GO:0031072", "GO:0043123", "GO:1901224", "GO:0043065",
        "GO:0051092", "GO:0007249", "GO:0090090", "GO:00431222", "GO:1902042")
go <- c("GO:2001240", "GO:0043065", "GO:2001238", "GO:2001234", "GO:0043524", "GO:2001243", "GO:0042771", "GO:0097193", "GO:1902042", "GO:0008625", "GO:0008625", 
        "GO:0006919", "GO:0051402", "GO:0006915", "GO:0007250", "GO:1901224", "GO:0043122","GO:0051059", "GO:0043123", "GO:0051092", "GO:0032088", "GO:0007249")
up <- RA_selected$Up_regulated
down <- RA_selected$Down_regulated
temp7 <- ""
up1 <- up1 <- c("EP300", "TP53BP2", "FAF1", "RNF25", "MTDH", "CHUK", "F2RL1", "IRAK1", "S100A4", "S100A13", "DDX21", "BRD4", "ANKRD17", "WLS", "TRIM8", "TRIM52", "MTDH", "DHX36", "JMJD8",
                "CAV1", "CHUK", "CLU", "DHX9", "EP300", "FLOT2", "HSPA1A", "HSPA1B", "IRAK1", "PLCG2", "PRKCI", "RPS6KA5", 
                "CLOCK", "ZBTB7A", "RNF25", "TRIM8", "TRIM52", "MTDH",
                "BMP7", "CHUK", "DDIT3", "IRAK1", "USP7", "TAX1BP1", "SIRT1", "ITCH", "OTULIN", "DAB2IP",
                "APAF1", "BAD", "HSPD1", "PDCD2", "EIF2AK3", "SENP1", "RPS27L", "DIABLO",
                "BCL2", "ATN1", "LIG4",
                "APAF1", "RHOB", "BAD", "BCL2",
"CASP3", "CD14", "DAPK3", "EP300", 
"FOXO1", "LTA", "MDM2", "MAP3K10", "DPF2", "STK4",
"PHLDA2", "FXR1", "OGT", "IER3", "BCLAF1", "THOC1", "BBC3", "DIABLO", "RNF152",
"AKT1", "CHUK", "ROCK1", "ROCK2", "IRAK1BP1",
"AKT1", "BCL2", "HSPA1A", "HSPA1B",
"IRAK1", "ZFP91",
"APAF1", "APC", "RHOB", 
"ATF4", "BAD", "BCL6", "BMP2", "C1QBP", "CLU", "CTNNB1", "DAPK3", "DDX3X", 
"FOXO1", "HTT", "FOXA1", "DNAJA1", "HSPD1", "ITGB1", "GADD45B", "NF1", "PDCD2", "PPID", "STK4", "HRK", "ARHGEF7", "SQSTM1", "SLIT2", "BCLAF1", "DNM1L", "RBM5", "RRP1B", "SIRT1", "RYBP", "SCRIB", "PHLDA3", "PNMA3", "TNFRSF12A", "ING3", "DIABLO", "SPHK2", "SUDS3", "ZNF622", "JMY", "DAB2IP", "C3orf38",
"CAV1", "PAK2", "TNFRSF12A",
"CD14", "DDX3X", "EP300", "IRAK1", "PHB", "NFAT5",
"BCL2", "CSNK2A2", "MAZ", "RB1", "MAPK8IP2",
"BCL2", "LIG4", "MSH2", "NAIP", "PRKCI", "SET", "TOX3", "FOXQ1",
"BCL2", "DDX3X", "NDUFS3", "URI1", "NOC2L",
"EP300", "SIRT1", "PHLDA3", "RPS27L", "AEN",
"F2RL1", "HSPB1", "PLCG2", "SQSTM1", "SPHK2", "DAB2IP", "SUMO4",
"APAF1", "BAD", "DDX3X", "CUL3", "CUL2", "DIABLO",
"DDX3X", "SERPINE1",
"BCL2", "DAXX", "DDX3X", "DIABLO", "DAB2IP")


up1 <- unique(up1)

saveRDS(up1, "upregulated_genes.rds")

#CORRELATED GENES TNFRSF12A
matrix<- GetAssayData(sce1, slot = "counts")
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["TNFRSF12A",])
correlations1<-apply(matrix_mod,1,function(x){cor(gene,x)})
corr <- sort(correlations1)
corr1 <- corr[31533:31633]
names_tnf <- names(corr1)

#New
png("Feb_Plots/bax.25_cluster1vs_pval1.png")
EnhancedVolcano(top1,
                lab = top1$gene,
                x = 'avg_log2FC',
                y = 'p_val',
                title = 'WT vs KO ALL DE Genes',
                pCutoff = .01,
                FCcutoff = .02,
                pointSize = 2.0,
                labSize = 3.0,
                xlim = c(-1, 1),
                col=c('blue', 'green', 'black', 'red3', "yellow", "brown", "grey"),
                colAlpha = 1)
dev.off()

Idents(sce1) <- sce1$cluster
sce1$cluster <- factor(sce1$cluster, level=c(1,2,3,4,5,6,7))
VlnPlot(sce1, features = c("TNFRSF12A"), pt.size = .1)

#subset cluster 3
sce1 <- subset(x = sce1, idents = "7")
sampnames4 <- substr(sce1@meta.data$Sample, 0, 2) 
Idents(sce1) <- sampnames4
all.markers <- unique(c(ofinterest, apotosis1, mDA.mature.markers, mDA.precursor.markers, A10.markers, A9.markers, non.mDA.markers))#, #ofinterest)) ,apotosis))
markers1 <- FindAllMarkers(sce1, min.pct = 0.001, min.diff.pct = 0.001, logfc.threshold = .001, features=tnf)
VlnPlot(sce1, features=all.markers)
top1 <- markers1 %>% dplyr::filter(cluster=="WT") 
#top <- top14[top14$gene %in% all.markers, ]
#top <- top14 %>% dplyr::filter(gene=c(all.markers))
genes <- top1 %>% dplyr::filter(p_val < .01)
png("Feb_Plots/WT_cluster3vln_final_final.png")
VlnPlot(sce1, features = genes$gene)
dev.off()




#GO ontologies
top1 <- markers1 %>% dplyr::filter(cluster=="WT")  %>% dplyr::filter(avg_log2FC > .01)
# Sort genes and only take positive genes



#Feature map
tnf <- c("TNFSF9", "TNFSF12", "TNFSF13B")
FeaturePlot(sce1, features = tnf)
a <- DotPlot(object =sce1, features = tnf)
a$data %>%
  kbl() %>%
  kable_styling()
