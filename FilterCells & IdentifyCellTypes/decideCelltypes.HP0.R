setwd("D:/scRNA-Seq/HP0_outs/")

load("HP0.no_doublets.Redblood.Rdata")
HP0 <- CreateSeuratObject(counts = HP0@assays$RNA@counts, min.cells = 10, min.features = 200, 
                          meta.data = HP0@meta.data)
HP0[["percent.mt"]] <- PercentageFeatureSet(HP0, pattern = "^MT-")
VlnPlot(HP0, features = c("nCount_RNA", 'nFeature_RNA', "percent.mt"), ncol = 3, pt.size = 0.2)
HP0 <- CellCycleScoring(HP0, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
HP0 <- NormalizeData(HP0, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T)
HP0 <- FindVariableFeatures(HP0, selection.method = "mean.var.plot")
HP0 <- ScaleData(HP0, vars.to.regress = 'Phase')
HP0 <- RunPCA(HP0, features = VariableFeatures(HP0), npcs = 20, ndims.print = 5, nfeatures.print = 5)
HP0 <- RunTSNE(object = HP0, dims = 1:12,check_duplicates = FALSE)

write.csv(colnames(HP0), 'HP0.barcodes.csv', row.names = FALSE, quote = FALSE)

markers=read.delim("D:/scRNA-Seq/aaaaaliteratures/2014.Blood.Global transcriptome analyses of human and murine terminal/edgeR/topMarkersForEachStage.Blood.txt",header=T,sep="\t",stringsAsFactors = F)
selected_markers=data.frame()
#select expressed genes and prepare proportion of gene expressing
for(s in colnames(markers)){
  stage=s
  select_markers=intersect(markers[,stage],rownames(HP0@assays$RNA@data))
  for(i in 1:length(select_markers)){
    tmp_markers=data.frame(stage=stage,Gene=select_markers[i],count=sum(HP0@assays$RNA@data[select_markers[i],]>0),rate=sum(HP0@assays$RNA@data[select_markers[i],]>0)/dim(HP0@assays$RNA@data)[2])
    selected_markers=rbind(selected_markers,tmp_markers)
  }
}
write.table(selected_markers,"topMarkersForEachStage.Blood.selected_markers.txt",sep="\t",row.names = F)

sigGenes.proerythroblast=selected_markers[selected_markers$stage=="proE"  & selected_markers$rate>0,"Gene"]
sigGenes.basophilic=selected_markers[selected_markers$stage=="basoE"  & selected_markers$rate>0,"Gene"]
sigGenes.polychromatic=selected_markers[selected_markers$stage=="polyE"  & selected_markers$rate>0,"Gene"]
sigGenes.orthochromatic=selected_markers[selected_markers$stage=="orthoE"  & selected_markers$rate>0,"Gene"]
modulegenes<- list(intersect(rownames(HP0@assays$RNA@data), sigGenes.proerythroblast))
HP0 <- AddModuleScore(HP0, features = modulegenes, name = "HP1.proerythroblast")
modulegenes<- list(intersect(rownames(HP0@assays$RNA@data), sigGenes.basophilic))
HP0 <- AddModuleScore(HP0, features = modulegenes, name = "HP1.basophilic")
modulegenes<- list(intersect(rownames(HP0@assays$RNA@data), sigGenes.polychromatic))
HP0 <- AddModuleScore(HP0, features = modulegenes, name = "HP1.polychromatic")
modulegenes<- list(intersect(rownames(HP0@assays$RNA@data), sigGenes.orthochromatic))
HP0 <- AddModuleScore(HP0, features = modulegenes, name = "HP1.orthochromatic")

HP0 <- RunTSNE(object = HP0, dims = 1:10,check_duplicates = FALSE)
FeaturePlot(HP0, cols = c("yellow", "red"), reduction = "tsne"
            , features =c("HP1.proerythroblast1","HP1.basophilic1","HP1.polychromatic1",'HP1.orthochromatic1')
            , pt.size = 0.5, ncol = 2)
FeaturePlot(HP0, cols = c("yellow", "red"), reduction = "tsne",
            features =c('GYPA', 'TFRC', 'HBG2','CENPA', 'TMCC2','CENPF'), pt.size = 0.5,ncol = 2)


# Determine cell stage
HP0 <- FindNeighbors(HP0, reduction = 'pca', dims = 1:10)
HP0 <- FindClusters(HP0,reduction.type="pca",dims.use=1:10,resolution = 0.6, print.output=0,save.SNN=T,force.recalc = T)
TSNEPlot(HP0, pt.size=0.2, label =TRUE)

uc_id<- as.numeric(as.vector(HP0@meta.data$RNA_snn_res.0.6))
uc_id[uc_id %in% c(8)]<-'proE'
uc_id[uc_id %in% c(6)]<-'basoE'
uc_id[uc_id %in% c(4)]<-'polyE'
uc_id[uc_id %in% c(0,1,2,3,5,7,9)]<-'orthoE'

clust_info<- data.frame(putative_stage=uc_id, row.names = rownames(HP0@meta.data), stringsAsFactors = FALSE)
HP0 <- AddMetaData(HP0, metadata = clust_info)

Idents(HP0) <- 'putative_stage'
HP0$putative_stage <- factor(HP0@meta.data$putative_stage, levels = c('proE','basoE','polyE','orthoE'))
TSNEPlot(HP0, group.by = 'putative_stage',pt.size=0.2,label = FALSE)
g = TSNEPlot(HP0, group.by = 'putative_stage',pt.size=0.2,label = FALSE)
table(ggplot_build(g)$data[[1]]$colour) # get the default colours
colors <- c( '#C77CFF','#F8766D','#00BFC4','#7CAE00')
TSNEPlot(HP0, group.by = 'putative_stage',pt.size=0.2) + 
  scale_color_manual(values = colors)

save(HP0,file="HP0.no_doublets.Redblood.stage.Rdata")

### evenly arrange Negative cells to two samples
load('D:/scRNA-Seq/HP0_outs/HP0.no_doublets.Redblood.stage.Rdata')
HP0@meta.data$HTO_classification <- as.character(HP0@meta.data$HTO_classification)
HP0@meta.data$HTO_classification_1 <- HP0@meta.data$HTO_classification
Negative_cells <- rownames(HP0@meta.data[HP0@meta.data$HTO_classification_1 == 'Negative',])
set.seed(1111)
hash1_cells <- sample(Negative_cells, round(length(Negative_cells)/2))
hash2_cells <- setdiff(Negative_cells, hash1_cells)
HP0@meta.data[hash1_cells, 'HTO_classification_1'] <- 'UCB2'
HP0@meta.data[hash2_cells, 'HTO_classification_1'] <- 'UCB3'
HP0@meta.data[HP0@meta.data$HTO_classification_1 == 'B0251 anti-human Hashtag1', 'HTO_classification_1'] <- 'UCB2'
HP0@meta.data[HP0@meta.data$HTO_classification_1 == 'B0252 anti-human Hashtag2', 'HTO_classification_1'] <- 'UCB3' 
save(HP0,file="HP0.no_doublets.Redblood.stage.Rdata")

colors <- c( '#C77CFF','#F8766D','#00BFC4','#7CAE00')
library(ggplot2)
TSNEPlot(HP0, group.by = 'putative_stage', pt.size = 0.5) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_color_manual(values = colors)
TSNEPlot(HP0, group.by = 'HTO_classification_1', pt.size = 0.5) + 
  theme_bw() +
  theme(panel.grid=element_blank())
TSNEPlot(HP0, group.by = 'Phase', pt.size = 0.5) + 
  theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c(G1='limegreen', G2M='blue', S= 'red'))

Idents(HP0) <- 'HTO_classification_1'
UCB2 <- subset(HP0, idents = 'UCB2')
UCB3 <- subset(HP0, idents = 'UCB3')
save(UCB2, file = 'UCB2.rdata')
save(UCB3, file = 'UCB3.rdata')

write.csv(colnames(UCB2), 'UCB2.barcodes.csv', row.names = FALSE, quote = FALSE)
write.csv(colnames(UCB3), 'UCB3.barcodes.csv', row.names = FALSE, quote = FALSE)


# re-clutering 
HP0 <- FindNeighbors(HP0, reduction = 'pca', dims = 1:10)
HP0 <- FindClusters(HP0,reduction.type="pca",dims.use=1:10,resolution = 0.1, print.output=0,save.SNN=T,force.recalc = T)
save(HP0,file="HP0.no_doublets.Redblood.stage.Rdata")
TSNEPlot(HP0, pt.size=0.2, label =TRUE)
markers <- FindAllMarkers(HP0, test.use = 'roc')
write.csv(markers, 'markers_for_5clusters.csv', row.names = FALSE, quote = FALSE)
