library(Rmisc)

load("HP1/h1_seurat.rdata")
RedbloodCells=rownames(seurat@meta.data[seurat@meta.data$Pred_lab=="Redblood",])
seurat=SubsetData(seurat,cells.use=RedbloodCells,subset.raw=T)

markers=read.delim("2014.Blood.Global transcriptome analyses of human and murine terminal/edgeR/topMarkersForEachStage.Blood.txt",header=T,sep="\t",stringsAsFactors = F)
selected_markers=data.frame()
#select expressed genes and prepare proportion of gene expressing
for(s in colnames(markers)){
  stage=s
  select_markers=intersect(markers[,stage],rownames(seurat@data))
  for(i in 1:length(select_markers)){
    tmp_markers=data.frame(stage=stage,Gene=select_markers[i],count=sum(seurat@data[select_markers[i],]>0),rate=sum(seurat@data[select_markers[i],]>0)/dim(seurat@data)[2])
    selected_markers=rbind(selected_markers,tmp_markers)
  }
}
write.table(selected_markers,paste0("HP1/HP1.","topMarkersForEachStage.Blood.selected_markers.txt"),sep="\t",row.names = F)

sigGenes.proerythroblast=selected_markers[selected_markers$stage=="proE"  & selected_markers$rate>0,"Gene"]
sigGenes.basophilic=selected_markers[selected_markers$stage=="basoE"  & selected_markers$rate>0,"Gene"]
sigGenes.polychromatic=selected_markers[selected_markers$stage=="polyE"  & selected_markers$rate>0,"Gene"]
sigGenes.orthochromatic=selected_markers[selected_markers$stage=="orthoE"  & selected_markers$rate>0,"Gene"]
modulegenes<- sigGenes.proerythroblast
modulegenes<- list(intersect(rownames(seurat@raw.data), modulegenes))
seurat <- AddModuleScore(seurat, genes.list = modulegenes, enrich.name = "HP1.proerythroblast")
modulegenes<- sigGenes.basophilic
modulegenes<- list(intersect(rownames(seurat@raw.data), modulegenes))
seurat <- AddModuleScore(seurat, genes.list = modulegenes, enrich.name = "HP1.basophilic")
modulegenes<- sigGenes.polychromatic
modulegenes<- list(intersect(rownames(seurat@raw.data), modulegenes))
seurat <- AddModuleScore(seurat, genes.list = modulegenes, enrich.name = "HP1.polychromatic")
modulegenes<- sigGenes.orthochromatic
modulegenes<- list(intersect(rownames(seurat@raw.data), modulegenes))
seurat <- AddModuleScore(seurat, genes.list = modulegenes, enrich.name = "HP1.orthochromatic")
FeaturePlot(seurat, cols.use = c("yellow", "red"), reduction.use = "tsne"
            , features.plot =c("HP1.proerythroblast1","HP1.basophilic1","HP1.polychromatic1",'HP1.orthochromatic1')
            , pt.size = 0.5, no.legend = F, max.cutoff = NA, dark.theme = T,nCol = 2)

#define cells of each stage
seurat=FindClusters(seurat,reduction.type="pca",dims.use=1:12,resolution=1.2,print.output=0,force.recalc=T)
TSNEPlot(seurat,do.label = T)
uc_id<- as.vector(seurat@ident)
plb<-uc_id
plb[uc_id %in% c(4)]<-'polyE'
plb[uc_id %in% c(0,1,2,3,5,6)]<-'orthoE'
clust_info<- data.frame(putative_stage=plb, row.names = seurat@cell.names,stringsAsFactors = F)
seurat<- AddMetaData(seurat, metadata =clust_info,col.name = "putative_stage")

basoE_cells=FeaturePlot(seurat,cols.use=c("yellow","red"),reduction.use="tsne",features.plot="HP1.basophilic1",do.identify = T,dark.theme =T)
seurat@meta.data[rownames(seurat@meta.data) %in% basoE_cells,"putative_stage"]<-"basoE" #58 cells

proE_cells=FeaturePlot(seurat,cols.use=c("yellow","red"),reduction.use="tsne",features.plot="HP1.proerythroblast1",do.identify = T,dark.theme =T)
seurat@meta.data[rownames(seurat@meta.data) %in% proE_cells,"putative_stage"]<-"proE" #5 cells

table(seurat@meta.data$putative_stage)
sum(table(seurat@meta.data$putative_stage))
dim(seurat@meta.data)
HP1_cell_summary=table(seurat@meta.data$putative_stage)
TSNEPlot(seurat,group.by="putative_stage",do.label=F,label.size = 8,no.legend=F)
save(seurat,file="HP1/h1_seurat.rdata")
