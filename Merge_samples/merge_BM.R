setwd("D:/scRNA-Seq/linwei/ModuleGenes/HuangPeng/")
load("HP2/h2_seurat.rdata")
HP2_seurat=seurat
load("HP3/h3_seurat.rdata")
HP3_seurat=seurat
load("HP4/h4_seurat.rdata")
HP4_seurat=seurat
HP2_seurat@meta.data$orig.ident="HP2"
HP2_seurat=RenameCells(HP2_seurat,add.cell.id="HP2")
HP3_seurat@meta.data$orig.ident="HP3"
HP3_seurat=RenameCells(HP3_seurat,add.cell.id="HP3")
HP4_seurat@meta.data$orig.ident="HP4"
HP4_seurat=RenameCells(HP4_seurat,add.cell.id="HP4")

HP23_seurat=MergeSeurat(HP2_seurat,HP3_seurat)
HP234_seurat=MergeSeurat(HP23_seurat,HP4_seurat)
rm(HP2_seurat,HP3_seurat,HP4_seurat,HP23_seurat)
cells.use=rownames(HP234_seurat@meta.data[is.na(HP234_seurat@meta.data$putative_stage)==F,])
HP234_seurat=SubsetData(HP234_seurat,cells.use=cells.use)

mito.features=grep(pattern="^MT-",rownames(HP234_seurat@raw.data), value = TRUE)
percent.mito=Matrix::colSums(HP234_seurat@raw.data[mito.features,])/Matrix::colSums(HP234_seurat@raw.data)
HP234_seurat=AddMetaData(HP234_seurat,metadata = percent.mito,col.name = "percent.mito")
VlnPlot(HP234_seurat,features.plot=c("percent.mito","nGene","nUMI"),nCol=3,group.by="putative_stage",do.sort=T,point.size.use=0.3)
VlnPlot(HP234_seurat,features.plot=c("percent.mito","nGene","nUMI"),nCol=3,group.by="orig.ident",do.sort=T,point.size.use=0.3)

HP234_seurat=NormalizeData(HP234_seurat,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP234_seurat=FindVariableGenes(HP234_seurat,mean.function=ExpMean,dispersion.function=LogVMR,
                        x.low.cutoff=0.0125,x.high.cutoff=3,y.cutoff=0.5,do.plot=F)
HP234_seurat<-ScaleData(HP234_seurat)
HP234_seurat<-RunPCA(HP234_seurat,pc.genes=HP234_seurat@var.genes,do.print=T,pcs.print=0,genes.print=0)
HP234_seurat=RunTSNE(HP234_seurat,do.fast=TRUE,dims.use=1:12)
save(HP234_seurat,file="HP234/HP234_seurat.Rdata")

HP234_seurat=SetAllIdent(HP234_seurat,id="putative_stage")
TSNEPlot(HP234_seurat,do.label=T,label.size=5,group.by="putative_stage")
TSNEPlot(HP234_seurat,do.label=T,label.size=5,group.by="orig.ident")
FeaturePlot(HP234_seurat,features.plot = "GYPA",no.legend = F)
FeaturePlot(HP234_seurat,features.plot = "MALAT1",no.legend = F)
FeaturePlot(HP234_seurat,features.plot = "NEAT1",no.legend = F)

HP234_seurat=FindClusters(HP234_seurat,dims.use = 1:12,force.recalc = T,resolution=0.5)
TSNEPlot(HP234_seurat,do.label=T,label.size=5)

table(HP234_seurat@meta.data$putative_stage)
table(HP234_seurat@meta.data$orig.ident)
HP234_seurat=FindClusters(HP234_seurat,dims.use = 1:12,force.recalc = T,resolution=3)
TSNEPlot(HP234_seurat,do.label=T,label.size=5) #for more clusters

HP234_seurat=CellCycleScoring(HP234_seurat, g2m.features=cc.genes$g2m.genes, s.features=cc.genes$s.genes, set.ident = FALSE)
TSNEPlot(HP234_seurat, group.by="Phase")

meta_data=HP234_seurat@meta.data
meta_data[meta_data$res.0.5 %in% c(4),"new_cluster"]="4"
meta_data[meta_data$res.0.5 %in% c(5),"new_cluster"]="5"
meta_data[meta_data$res.0.5 %in% c(1,2,6),"new_cluster"]="126"
meta_data[meta_data$res.0.5 %in% c(3),"new_cluster"]="3"
meta_data[meta_data$res.0.5 %in% c(0),"new_cluster"]="0"
meta_data[meta_data$res.0.5 %in% c(7),"new_cluster"]="7"
HP234_seurat@meta.data=meta_data
HP234_seurat=SetAllIdent(HP234_seurat,id="new_cluster")
markers_cluster=FindAllMarkers(HP234_seurat,only.pos=T,test.use="roc")
write.table(markers_cluster,"HP234/HP234_seurat.markers_cluster.roc.txt",sep="\t",row.names=F)

library(dplyr)
top_markers_cluster=markers_cluster %>% group_by(cluster) %>% top_n(30,-p_val)
top_markers_cluster$gene=as.character(top_markers_cluster$gene)
genes_plot=c(top_markers_cluster$gene[18:284],top_markers_cluster$gene[315:344],
             top_markers_cluster$gene[1:17],top_markers_cluster$gene[285:314])
#HBG1:17,TMSB4X:285,TYROBP:314
genes_plot=c(genes_plot[1:181],genes_plot[208:267],genes_plot[182:207],genes_plot[268:344])
DoHeatmap(HP234_seurat,genes.use=genes_plot,group.order = c(3,4,5,126,0,7))

#GO analysis
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
go=read.gmt(gmtfile="E:/scRNA-Seq/global_data/refdata/c5.all.v6.2.symbols.gmt")
genes_used=markers_cluster[markers_cluster$p_val<0.05 & markers_cluster$cluster=="7","gene"]
eg=bitr(genes_used, fromType="SYMBOL", toType=c("SYMBOL","ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
eg_go=enricher(eg$SYMBOL,TERM2GENE=go,pvalueCutoff = 0.1)
dotplot(eg_go,title='GO',showCategory=30,font.size=8)
write.table(summary(eg_go),"HP234/eg_go.cluster7.txt",row.names = F,sep="\t")

#look at cluster 4,5
c4c5_seurat=SubsetData(HP234_seurat,ident.use=c(4,5))
markers_cluster_c4c5=FindAllMarkers(c4c5_seurat,only.pos=T,test.use="roc")
write.table(markers_cluster_c4c5,"HP234/HP234_seurat.markers_cluster_c4c5.roc.txt",sep="\t",row.names=F)


# throw cluster 7 away
Idents(HP234_seurat) <- 'new_cluster'
HP234_seurat <- subset(HP234_seurat, idents = 7 ,invert = TRUE)
save(HP234_seurat, file = 'HP234/HP234_seurat.Rdata')
