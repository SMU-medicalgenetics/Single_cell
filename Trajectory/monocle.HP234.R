library(monocle)
library(Seurat)
load("HP234/HP234_seurat.Rdata")

HSMM=importCDS(HP234_seurat,import_all=T)
HSMM=estimateSizeFactors(HSMM)
HSMM=estimateDispersions(HSMM)
HSMM=detectGenes(HSMM,min_expr=0)

#Clustering cells without marker genes
disp_table=dispersionTable(HSMM)
unsup_clustering_genes=subset(disp_table,mean_expression>=0) #no filtering

#Constructing Single Cell Trajectories
HSMM=setOrderingFilter(HSMM,unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
HSMM=reduceDimension(HSMM,max_components=2,num_dim=20,reduction_method='DDRTree',verbose=F,norm_method = "none")
HSMM=orderCells(HSMM,reverse=T)
pData(HSMM)$putative_stage=HP234_seurat@meta.data[rownames(pData(HSMM)),"putative_stage"]
save(HSMM,file="HP234_HSMM.Rdata")
plot_cell_trajectory(HSMM,color_by="Pseudotime",markers = "")
plot_cell_trajectory(HSMM,color_by="putative_stage",markers = "")
colnames(pData(HSMM))[63]="Cluster" #res.0.5 to Clusetr
plot_cell_trajectory(HSMM,color_by="Cluster",markers = "")

hm=plot_pseudotime_heatmap(HSMM,num_clusters=4,show_rownames=F,return_heatmap = T)
save(hm,file="HP234.heatmap.Rdata")
c1 <- as.data.frame(cutree(hm$tree_row, k=4)) 
colnames(c1) <- "Cluster"
c1$Gene <- rownames(c1)
write.table(c1,"HP234.genes_in_heatmap_clusters.txt",sep="\t",row.names = F)

#plot some genes
genes_1=c("CFL1","SVIP","EIF5","SLC2A1")
genes_2=c("RNF10", "TERF2IP", "AKAP8L")
genes_3=c("POLR2L", "YBX1")

library(smoother)
library(plotrix)
library(ggpubr)
plots <- list()
info=pData(HSMM)
info=info[order(info$Pseudotime),]
ggplot(info,aes(x=`TSNE-1`,y=`TSNE-2`,color=Pseudotime))+geom_point()
x_axis=1:dim(info)[1]
for(gene_id in genes_3){
  y_axis <- as.numeric(exprs(HSMM[gene_id,rownames(info)]))
  lo <- loess(y_axis~x_axis)
  xl <- seq(min(x_axis),max(x_axis),(max(x_axis) - min(x_axis))/1000)
  y_axis <- predict(lo,xl)
  #y_axis <- rescale(y_axis, newrange = c(0,1))
  df <- data.frame(cells=1:length(y_axis),Expression=as.numeric(y_axis))
  df$gene <- gene_id
  df$cells <- factor(df$cells, levels=df$cells)
  num_subpops <- length(unique(df$population))
  plots[[gene_id]] <- df
}
tables <- do.call(rbind, plots)
tables$trajectory="o"
alldfs <- tables

ggplot(alldfs, aes(x=cells, y=Expression, group=gene, colour=gene)) +
  theme_bw() + geom_line(size=1) + xlab('Pseudotime') +
  guides(col=guide_legend(direction="vertical")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks=element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5))+
  scale_color_manual("", values=c("purple","gold3","cyan","blue","red","darkgreen","yellow","black","chocolate1")) 
