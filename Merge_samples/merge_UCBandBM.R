setwd('D:/scRNA-Seq/HP01234/')
load('D:/scRNA-Seq/linwei/ModuleGenes/HuangPeng/HP234/HP234_seurat.Rdata')
load('D:/scRNA-Seq/HP0_outs/HP0.no_doublets.Redblood.stage.Rdata')
HP0@meta.data$orig.ident <- HP0@meta.data$HTO_classification_1
CBBM <- merge(HP0, HP234_seurat, merge.data = TRUE)
CBBM@meta.data[CBBM@meta.data$orig.ident %in% c('HP2','HP3','HP4'), 'source'] <- 'BM'
CBBM@meta.data[CBBM@meta.data$orig.ident %in% c('UCB2','UCB3'), 'source'] <- 'UCB'
  
CBBM <- NormalizeData(CBBM, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
CBBM <- ScaleData(CBBM, vars.to.regress = 'orig.ident')
CBBM <- FindVariableFeatures(CBBM, selection.method = "vst", nfeatures = 800)
CBBM <- RunPCA(CBBM, features = VariableFeatures(CBBM))
CBBM <- RunTSNE(object = CBBM, dims = 1:20,check_duplicates = FALSE)
save(CBBM, file = 'CBBM.rdata')

library(ggplot2)
TSNEPlot(CBBM, group.by = 'source') + theme_bw() +
  theme(panel.grid=element_blank())
TSNEPlot(CBBM, group.by = 'orig.ident') + theme_bw() +
  theme(panel.grid=element_blank())
TSNEPlot(CBBM, group.by = 'putative_stage') + theme_bw() +
  theme(panel.grid=element_blank())
FeaturePlot(CBBM, features = hbgenes, slot = 'data')


### plot Doheatmap and correlations between samples 
Idents(CBBM) <- 'source'
markers <- FindAllMarkers(CBBM,test.use = 'wilcox', only.pos = T)
topgene <- markers %>% group_by(cluster) %>% top_n(200,  avg_logFC)
DoHeatmap(CBBM, topgene$gene, size=3, group.by = 'orig.ident')

# explore why separate UCB and BM, even throw hb genes and try again
markers <- FindAllMarkers(CBBM,test.use = 'wilcox', only.pos = T, min.pct = 0.1, min.diff.pct = 0.25)
write.csv(markers, 'markers.source.csv', row.names = F)
FeaturePlot(CBBM, features = markers[markers$cluster == 'BM', 'gene'], slot = 'data')
FeaturePlot(CBBM, features = markers[markers$cluster == 'UCB', 'gene'][1:9], slot = 'data')

hbgenes <- c('HBG2','HBG1',"HBB","HBD","HBZ","HBA1","HBA2","HBM","HBP1","HBQ1")
keep.genes <- setdiff(rownames(CBBM@assays$RNA@counts), hbgenes)
CBBM_sub <- CBBM[keep.genes, ]
CBBM_sub <- NormalizeData(CBBM_sub, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
CBBM_sub <- ScaleData(CBBM_sub, vars.to.regress = 'orig.ident')
CBBM_sub <- FindVariableFeatures(CBBM_sub, selection.method = "vst", nfeatures = 800)
CBBM_sub <- RunPCA(CBBM_sub, features = VariableFeatures(CBBM_sub))
CBBM_sub <- RunTSNE(object = CBBM_sub, dims = 1:20,check_duplicates = FALSE)
save(CBBM_sub, file = 'CBBM_sub.rdata')

TSNEPlot(CBBM_sub, group.by = 'source') + theme_bw() +
  theme(panel.grid=element_blank())
TSNEPlot(CBBM_sub, group.by = 'orig.ident') + theme_bw() +
  theme(panel.grid=element_blank())
TSNEPlot(CBBM_sub, group.by = 'putative_stage') + theme_bw() +
  theme(panel.grid=element_blank())
FeaturePlot(CBBM_sub, features = hbgenes, slot = 'data')

VizDimLoadings(CBBM, dims = 1:10)
DimHeatmap(CBBM, dims = 1:10, cells = 500, balanced = TRUE, ncol = 3)