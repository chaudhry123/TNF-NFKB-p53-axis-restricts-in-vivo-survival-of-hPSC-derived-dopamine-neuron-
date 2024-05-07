#!/usr/bin/env Rscript

# set up input directory
inputd = "TNF-NFKB-p53-axis-restricts-in-vivo-survival-of-hPSC-derived-dopamine-neuron-/Bulk-RNAseq-processing/internal-data"

library(data.table)
library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(NMF)
library(amap)
library(dendextend)
library(msigdbr)
library(dplyr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(enrichplot)
library(openxlsx)
library(ggpubr)
library(magrittr)

fn = dir(path=inputd, pattern="*-htseq-count-matrix.txt", recursive=T, full.names=T)
fn %>% file.exists
# [1] TRUE TRUE TRUE TRUE TRUE TRUE

mySampleName = gsub("Proj_11351_C_s_E223", "", gsub("-htseq.*", "", basename(fn)))
mySampleTable = data.frame(sampleName=mySampleName, fileName=fn, condition=factor(c("D0", "D0", "D1Culture", "D1Culture", "D1Graft", "D1Graft"), levels=c("D0", "D1Culture", "D1Graft")))

DESeq2Obj = DESeqDataSetFromHTSeqCount(sampleTable=mySampleTable, directory=".", design= ~ condition)
DESeq2Obj %<>% DESeq
DESeq2Obj.vst = varianceStabilizingTransformation(DESeq2Obj)
dim(DESeq2Obj)
# [1] 20492     6

p1 = plotPCA(DESeq2Obj.vst, ntop=3000, intgroup="condition", returnData=T)

percentVar = round(100 * attr(p1, "percentVar"))
p1.4a = ggplot() + geom_point(data=p1, aes(PC1, PC2, color=condition), size=4.0) +
  xlab(paste0("PC1: ", percentVar[1], "%")) + ylab(paste0("PC2: ", percentVar[2], "%")) + coord_fixed() +
  theme(panel.grid=element_blank()) +
  geom_text_repel(data=p1, aes(PC1, PC2, label=row.names(p1)), size=3, max.overlaps=30, segment.alpha=0.1, box.padding=0.2, point.padding=0.20)

p1
#                        PC1           PC2     group condition            name
# D0P531WT1        -5.703104  55.192211708        D0        D0       D0P531WT1
# D0P531WT2        37.338385  30.414625140        D0        D0       D0P531WT2
# D1CULTUREP53WT1  33.434946 -23.298629845 D1Culture D1Culture D1CULTUREP53WT1
# D1CULTUREP53WT2  22.213052 -14.094382774 D1Culture D1Culture D1CULTUREP53WT2
# D1GRAFTP53WT1    -5.617410 -48.209027546   D1Graft   D1Graft   D1GRAFTP53WT1
# D1GRAFTP53WT2   -81.665868  -0.004796684   D1Graft   D1Graft   D1GRAFTP53WT2

ntop.3000.selector = rowVars(assay(DESeq2Obj.vst)) %>% order(decreasing=T) %>% head(3000)
row.hclust = hcluster(dist(t(assay(DESeq2Obj.vst)[ntop.3000.selector, ])), method="euclidean", link="single")
row.hclust.d = as.dendrogram(row.hclust)
labels(row.hclust.d)
# [1] "D1GRAFTP53WT2"   "D1GRAFTP53WT1"   "D1CULTUREP53WT1" "D1CULTUREP53WT2"
# [5] "D0P531WT1"       "D0P531WT2"

# normalized DESeq2 counts, log2 transform, top 3000 genes, correlation metric, ward linkage
row.hclust.all = hcluster(dist(t(log2(counts(DESeq2Obj, normalized=T) + 1)[ntop.3000.selector, ])), method="correlation", link="ward")
row.hclust.all.d = as.dendrogram(row.hclust.all)
row.hclust.all.d = row.hclust.all.d %>% set("labels", c("Day 1 Culture replicate 1", "Day 1 Culture replicate 2", "Day 0 replicate 1", "Day 0 replicate 2", "Day 1 Graft replicate 1", "Day 1 Graft replicate 2")) %>% set("leaves_pch", 19) %>% set("leaves_cex", 1.5) %>% set("leaves_col", c("green", "green", "red", "red", "blue", "blue"))
par(mar=c(0.6, 0.6, 0.6, 10.1))
pdf("figures4a.pdf", width=7, height=7)
print(p1.4a)
plot(row.hclust.all.d, sub="", xlab="", yaxt="n", horiz=T)
dev.off()

diffExprLst.list = list(
  D0.D1Culture=results(DESeq2Obj, contrast=c("condition", "D1Culture", "D0"), alpha=0.05),
  D0.D1Graft=results(DESeq2Obj, contrast=c("condition", "D1Graft", "D0"), alpha=0.05),
  D1Culture.D1Graft=results(DESeq2Obj, contrast=c("condition", "D1Graft", "D1Culture"), alpha=0.05)
)
for (j in 1:length(diffExprLst.list)) {
  diffExprLst.list[[j]]$padj[is.na(diffExprLst.list[[j]]$padj)] = 1.0
}

p1 = plotMA(diffExprLst.list[[3]], alpha=0.05, returnData=T)
pdf("figure4b.pdf", width=7, height=7)
plotMA(diffExprLst.list[[3]], alpha=0.05, main="D1 graft vs D1 culture", sub="Q value 0.05", xlab="", ylim=c(-6.0, +6.0), ylab="log2-fold change", cex.axis=0.80, colSig=NA)
points(pmin(lfc, +6.0) ~ mean, data=subset(p1, isDE & lfc > 0), col="red", pch=ifelse(lfc >= +6.0, 2, 20), cex=0.45)
points(pmax(lfc, -6.0) ~ mean, data=subset(p1, isDE & lfc < 0), col="blue", pch=ifelse(lfc <= -6.0, 6, 20), cex=0.45)
title(xlab="mean of normalized counts", cex.lab=0.80)
legend("topright", inset=c(0, 0), legend=sum(diffExprLst.list[[3]]$padj <= 0.05 & diffExprLst.list[[3]]$log2FoldChange >= 0), bty="n", cex=1.30, text.col="red")
legend("bottomright", inset=c(0, 0), legend=sum(diffExprLst.list[[3]]$padj <= 0.05 & diffExprLst.list[[3]]$log2FoldChange < 0), bty="n", cex=1.30, text.col="blue")
dev.off()

packageVersion("msigdbr")
# [1] '7.5.1'
h_gene_sets = msigdbr(species="human", category="H")
h_gene_sets$gs_name = gsub("^HALLMARK_", "", h_gene_sets$gs_name)

msigdbr_t2g = h_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
hallmark.df.3.up = enricher(gene=rownames(subset(diffExprLst.list[[3]], log2FoldChange >= 0 & padj <= 0.05)), TERM2GENE=msigdbr_t2g)
hallmark.df.3.down = enricher(gene=rownames(subset(diffExprLst.list[[3]], log2FoldChange < 0 & padj <= 0.05)), TERM2GENE=msigdbr_t2g)

pdf("figure4c.pdf", width=7, height=7)
dotplot(hallmark.df.3.up, showCategory=30) + theme_classic()
dev.off()
pdf("figures4b.pdf", width=7, height=7)
dotplot(hallmark.df.3.down, showCategory=30) + theme_classic()
dev.off()

go.term = "GO:0034612" # (response to tumor necrosis factor)

# The publication used version 3.13 of org.Hs.eg.db.
packageVersion("org.Hs.eg.db")
# [1] '3.19.1'

go_allegs_object = as.list(org.Hs.egGO2ALLEGS)
# Difference between org.Hs.egGO2EG and org.Hs.eg.GO2ALLEGS
# ALLEGS finds genes that are associated up to child nodes.
entrez.fig.s5c = go_allegs_object$`GO:0034612` %>% unique
symbols.figs.s5c = AnnotationDbi::select(org.Hs.eg.db, keys=entrez.fig.s5c, columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID")
length(symbols.figs.s5c$SYMBOL)


heatmap.input.withvar = assay(DESeq2Obj.vst)[rowVars(counts(DESeq2Obj)) > 0, ]
heatmap.gsea.withvar.1 = pheatmap:::scale_rows(heatmap.input.withvar[rownames(heatmap.input.withvar) %in% symbols.figs.s5c$SYMBOL, ] %>% as.matrix)

# Placed some cutoff on gene expression, on day 1 graft replicates 1 and 2.
gsea.56up.genes.list = rownames(heatmap.gsea.withvar.1[heatmap.gsea.withvar.1[, 5] >= 0.3 & heatmap.gsea.withvar.1[, 6] >= 0.3, ])

pheatmap(heatmap.input.withvar[rownames(heatmap.input.withvar) %in% (gsea.56up.genes.list %>% unique), ] %>% as.matrix, scale="row", main="GO:0034612 (response to tumor necrosis factor)", color=colorRampPalette(c("blue", "white", "red"))(100), cluster_rows=T, cluster_distance_rows="euclidean", cluster_cols=F, show_rownames=T, fontsize=7, filename="figs4c-2.pdf")

# GSEA using gseaplot2
gene.list.for.gsea = diffExprLst.list[[3]]$log2FoldChange
names(gene.list.for.gsea) = rownames(diffExprLst.list[[3]])
gene.list.for.gsea = na.omit(gene.list.for.gsea)
gene.list.for.gsea = sort(gene.list.for.gsea, decreasing=T)

gse = gseGO(gene.list.for.gsea, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP", pvalueCutoff=1.0)
# this part is not reproducible. Saving the result is recommended.
dim(gse)
# [1] 6153   11

# GSEA plot using gseaplot2
j = grep("tumor necro", gse$Description)[gse[grep("tumor necro", gse$Description), ]$p.adjust %>% which.min]
message(j)

gse[j, c("ID", "Description", "NES", "p.adjust")]
#                    ID                       Description      NES   p.adjust
# GO:0034612 GO:0034612 response to tumor necrosis factor 1.474993 0.06298552

pdf("figure4d.pdf", width=7, height=7)
gseaplot2(gse, title=gse$Description[j], geneSetID=j)
dev.off()

gse.ordered = read.xlsx("internal-data/GSEA.1c1g.ranked.by.NES.6181.terms.20211122_highlight.xlsx")
all.equal(gse.ordered$NES, gse.ordered$NES %>% sort(decreasing=T))
# [1] TRUE

# New criterion here: top 500. It is justified because the following edge (rank 1401 - 1500) has NES around 1.1 or 1.2, low.
gse.ordered.1 = gse.ordered[1:500, ]
dim(gse.ordered.1)
# [1] 500  11
all.equal(gse.ordered.1$NES, gse.ordered$NES[1:500])
# [1] TRUE
all.equal(gse.ordered.1$NES, gse.ordered.1$NES %>% sort(decreasing=T))
# [1] TRUE
nes.annotated = read.xlsx("GSEA rank plot  FINAL list.xlsx", sheet=1, colNames=F)
nes.annotated$X2 %in% gse.ordered.1$ID %>% all
# [1] TRUE
gse.ordered.1$color.selector = ifelse(gse.ordered.1$ID %in% nes.annotated$X2, "pass", "nopass")
gse.ordered.1$rankx = 1:NROW(gse.ordered.1)
gse.ordered.1$Legend = "nopass"
gse.ordered.1$Legend[gse.ordered.1$ID %in% nes.annotated[1:6, ]$X2] = "TNF alpha"
gse.ordered.1$Legend[gse.ordered.1$ID %in% nes.annotated[7:12, ]$X2] = "Apoptosis"
gse.ordered.1$Legend[gse.ordered.1$ID %in% nes.annotated[13:16, ]$X2] = "P53"
gse.ordered.1$pointsize = ifelse(gse.ordered.1$Legend == "nopass", 3, 8)

gse.ordered.1 = gse.ordered.1[order(gse.ordered.1$color.selector), ]
# Alphabetical order. The order is "nopass" then "pass", which is what we want.

gse.ordered.1$Legend %>% table
# Apoptosis    nopass       P53 TNF alpha 
#         6       484         4         6
# The highlighted GO terms, manually curated
highlighted.gsea = list(
  `TNF alpha`=nes.annotated[1:6, 2:3],
  `Apoptosis`=nes.annotated[7:12, 2:3],
  `P53`=nes.annotated[13:16, 2:3]
)
highlighted.gsea
# $`TNF alpha`
#           X2
# 1 GO:0034612
# 2 GO:0071356
# 3 GO:0032680
# 4 GO:0032640
# 5 GO:0071706
# 6 GO:1903555
#                                                                    X3
# 1                                   response to tumor necrosis factor
# 2                          cellular response to tumor necrosis factor
# 3                      regulation of tumor necrosis factor production
# 4                                    tumor necrosis factor production
# 5               tumor necrosis factor superfamily cytokine production
# 6 regulation of tumor necrosis factor superfamily cytokine production
# 
# $Apoptosis
#            X2
# 7  GO:0097192
# 8  GO:2001239
# 9  GO:1902229
# 10 GO:0097191
# 11 GO:1900117
# 12 GO:2001242
#                                                                               X3
# 7                     extrinsic apoptotic signaling pathway in absence of ligand
# 8       regulation of extrinsic apoptotic signaling pathway in absence of ligand
# 9  regulation of intrinsic apoptotic signaling pathway in response to DNA damage
# 10                                         extrinsic apoptotic signaling pathway
# 11                                    regulation of execution phase of apoptosis
# 12                           regulation of intrinsic apoptotic signaling pathway
# 
# $P53
#            X2
# 13 GO:0043516
# 14 GO:1902165
# 15 GO:1902253
# 16 GO:1901798
#                                                                                                     X3
# 13                        regulation of DNA damage response, signal transduction by p53 class mediator
# 14 regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator
# 15                           regulation of intrinsic apoptotic signaling pathway by p53 class mediator
# 16                                    positive regulation of signal transduction by p53 class mediator
# Add "GO: " before the term names
gse.ordered.1$Description = paste0("GO: ", gse.ordered.1$Description)

nes.1 = ggplot() + geom_point(data=gse.ordered.1, aes(rankx, abs(NES), fill=Legend, size=pointsize), shape=21, color=rgb(0, 0, 0, 0.10), show.legend=F) + scale_size_identity() + labs(x="Rank", y="Absolute value of GSEA NES") + scale_fill_manual(values=c("TNF alpha"="red", "Apoptosis"="saddlebrown", "P53"="purple", "nopass"=rgb(0, 0, 0, 0.10))) + theme_classic() + theme(plot.title=element_text(hjust=0.5)) + geom_text_repel(data=subset(gse.ordered.1, Legend != "nopass"), aes(rankx, abs(NES), label=subset(gse.ordered.1, Legend != "nopass")$Description, color=Legend), key_glyph="rect", size=4, max.overlaps=100, segment.alpha=0.2, box.padding=0.4, point.padding=1, force=2000, max.iter=3000000) + scale_color_manual(values=c("TNF alpha"="red", "Apoptosis"="saddlebrown", "P53"="purple"))

pdf("figure4e.pdf", width=7, height=7)
print(nes.1)
dev.off()


library(devtools)
session_info(pkgs="attached")

