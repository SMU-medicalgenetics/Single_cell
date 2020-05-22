setwd('D:/scRNA-Seq/HP0_outs/')
library(Seurat)
data_dir <- "D:/scRNA-Seq/HP0_outs/outs/filtered_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir)
HP0 <- CreateSeuratObject(counts = data$`Gene Expression`)
HP0 <- NormalizeData(HP0)
HP0 <- FindVariableFeatures(HP0, selection.method = "mean.var.plot")
HP0 <- ScaleData(HP0, features = VariableFeatures(HP0))
HP0[['HTO']] = CreateAssayObject(counts = data$`Antibody Capture`)
HP0 <- NormalizeData(HP0, assay = "HTO", normalization.method = "CLR")
#Demultiplex cells based on HTO enrichment to assign single cells back to their sample origins.
HP0 <- HTODemux(HP0, assay = "HTO", positive.quantile = 0.85, kfunc = 'kmeans')

Idents(HP0) <- "HTO_maxID"
RidgePlot(HP0, assay = "HTO", features = rownames(HP0[["HTO"]]), ncol = 2)
FeatureScatter(HP0, feature1 = "B0251 anti-human Hashtag1", feature2 = "B0252 anti-human Hashtag2")
Idents(HP0) <- "HTO_classification.global"
VlnPlot(HP0, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

# Strategy: all cells but doublets of two samples to predict cell types using SuperCT framework, only keep Redblood cells
Idents(HP0) <- "HTO_classification.global"
HP0 <- subset(HP0, idents = c('Negative','Singlet'))

library(rSuperCT)
library(ggplot2)
pred_obj <- ImportData(HP0)
pred_obj <- PredCellTypes(pred_obj, model = 'CellTypes', results.dir = '~/rSuperCT_Files/TestFiles')
g <- plotHist(pred_obj) + scale_fill_manual(values = rep('blue', 13))
write.csv(g$data, 'pred_types.hist.csv',row.names = F)

HP0@meta.data$pred_types <- pred_obj@meta.data$pred_types
table(HP0@meta.data$pred_types,HP0@meta.data$HTO_classification.global)

#20200407: split samples backward
# Caution: use classification of HP0 in HP0.no_doublets.Redblood.stage.Rdata
HP0@meta.data$HTO_classification_1 <- NULL
HP0@meta.data$HTO_classification_1 <- meta.data[rownames(HP0@meta.data), 'HTO_classification_1']
pred_types <- unique(HP0@meta.data$pred_types)
set.seed(111)
for(t in pred_types){
  tmp_cells <- rownames(HP0@meta.data[HP0@meta.data$pred_types == t & is.na(HP0@meta.data$HTO_classification_1),])
  cells_one <- sample(tmp_cells, round(length(tmp_cells)/2))
  cells_two <- setdiff(tmp_cells, cells_one)
  HP0@meta.data[cells_one, 'HTO_classification_1'] <- 'UCB2'
  HP0@meta.data[cells_two, 'HTO_classification_1'] <- 'UCB3'
}
table(HP0$pred_types,HP0$HTO_classification_1)

Idents(HP0) <- "pred_types"
HP0 <- subset(HP0, idents = 'Redblood')
save(HP0, file = 'HP0.no_doublets.Redblood.Rdata')

