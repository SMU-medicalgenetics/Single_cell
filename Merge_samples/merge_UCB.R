setwd('D:/scRNA-Seq/HP01')
load('D:/scRNA-Seq/HP1_outs/outs/HP1P_outs/HP1P.Rdata')
HP1 <- UpdateSeuratObject(HP1)
HP1@meta.data$orig.ident <- 'UCB1'
UCB1 <- HP1
colnames(UCB1@assays$RNA@counts) <- paste0(colnames(UCB1@assays$RNA@counts),"_1")
rownames(UCB1@meta.data) <- colnames(UCB1@assays$RNA@counts)
UCB1 <- CreateSeuratObject(counts = UCB1@assays$RNA@counts, meta.data = UCB1@meta.data)
UCB1 <- FindVariableFeatures(UCB1, selection.method = "vst", nfeatures = 2000)
UCB1 <- NormalizeData(UCB1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
UCB1 <- ScaleData(UCB1)
UCB1 <- RunPCA(UCB1, features = VariableFeatures(UCB1), npcs = 20)
UCB1 <- RunTSNE(object = UCB1, dims = 1:12,check_duplicates = FALSE)

UCB1 <- CellCycleScoring(UCB1, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
TSNEPlot(UCB1 , group.by='Phase')

phase <- read.csv('D:/scRNA-Seq/HP1_outs/scImpute_k5/cellrangerRkit/Graph-based.with_clusterID.with_stage.csv', stringsAsFactors = FALSE)
rownames(phase) <- gsub('-1', '_1', phase$Barcode)
UCB1@meta.data$Phase <- phase[rownames(UCB1@meta.data), 'phases']
TSNEPlot(UCB1 , group.by='Phase')

load('D:/scRNA-Seq/HP0_outs/UCB2.rdata')
UCB2@meta.data$orig.ident <- 'UCB2'
UCB2@assays$HTO <- NULL
colnames(UCB2@assays$RNA@counts) <- paste0(colnames(UCB2@assays$RNA@counts),"_2")
rownames(UCB2@meta.data) <- colnames(UCB2@assays$RNA@counts)
UCB2 <- CreateSeuratObject(counts = UCB2@assays$RNA@counts, meta.data = UCB2@meta.data)
UCB2 <- FindVariableFeatures(UCB2, selection.method = "vst", nfeatures = 2000)
UCB2 <- NormalizeData(UCB2, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
UCB2 <- ScaleData(UCB2)

load('D:/scRNA-Seq/HP0_outs/UCB3.rdata')
UCB3@meta.data$orig.ident <- 'UCB3'
UCB3@assays$HTO <- NULL
colnames(UCB3@assays$RNA@counts) <- paste0(colnames(UCB3@assays$RNA@counts),"_3")
rownames(UCB3@meta.data) <- colnames(UCB3@assays$RNA@counts)
UCB3 <- CreateSeuratObject(counts = UCB3@assays$RNA@counts, meta.data = UCB3@meta.data)
UCB3 <- FindVariableFeatures(UCB3, selection.method = "vst", nfeatures = 2000)
UCB3 <- NormalizeData(UCB3, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
UCB3 <- ScaleData(UCB3)

UCB12 <- RunCCA(object1=UCB1,object2=UCB2)
UCB12 <- FindVariableFeatures(UCB12, selection.method = "vst", nfeatures = 2000)
UCB123 <- RunCCA(object1=UCB12,object2=UCB3)

# UCB123@meta.data[rownames(UCB1@meta.data), 'Phase'] <- as.character(UCB1@meta.data$Phase)
# UCB123@meta.data[rownames(UCB2@meta.data), 'Phase'] <- as.character(UCB2@meta.data$Phase)
# UCB123@meta.data[rownames(UCB3@meta.data), 'Phase'] <- as.character(UCB3@meta.data$Phase)

UCB123 <- NormalizeData(UCB123, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
UCB123 <- ScaleData(UCB123, vars.to.regress = 'orig.ident')
UCB123 <- FindVariableFeatures(UCB123, selection.method = "vst", nfeatures = 2000)
UCB123 <- RunPCA(UCB123, features = VariableFeatures(UCB123), npcs = 20)
UCB123 <- RunTSNE(object = UCB123, dims = 1:12,check_duplicates = FALSE)

UCB123$putative_stage <- factor(UCB123@meta.data$putative_stage, levels = c('proE','basoE','polyE','orthoE'))
save(UCB123, file = 'UCB123.rdata')

colors <- c( '#C77CFF','#F8766D','#00BFC4','#7CAE00')
library(ggplot2)
TSNEPlot(UCB123, group.by = 'putative_stage', pt.size = 0.5) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_color_manual(values = colors)
TSNEPlot(UCB123, group.by = 'orig.ident', pt.size = 0.5) + 
  theme_bw() +
  theme(panel.grid=element_blank())
TSNEPlot(UCB123, group.by = 'Phase', pt.size = 0.5) + 
  theme_bw() +
  theme(panel.grid=element_blank())+
  scale_color_manual(values=c(G1='limegreen', G2M='blue', S= 'red'))


# basic stat
median(UCB123@meta.data$nFeature_RNA)
mean(UCB123@meta.data$nFeature_RNA)
aggregate(nFeature_RNA~orig.ident, UCB123@meta.data, median)

median(UCB123@meta.data$nCount_RNA)
mean(UCB123@meta.data$nCount_RNA)
aggregate(nCount_RNA~orig.ident, UCB123@meta.data, median)
