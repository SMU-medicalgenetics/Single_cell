# Seurat version 2.3.4
#use scImpute before downstream analysis
library(scImpute)
library(Seurat)
setwd("E:/scRNA-Seq/HP4_outs/scImpute_k5")
HP4.data<-Read10X(data.dir="E:/scRNA-Seq/HP4_outs/outs/filtered_gene_bc_matrices/hg19")
write.table(as.matrix(HP4.data),"E:/scRNA-Seq/HP4_outs/HP4_data.txt",sep="\t")
scimpute("E:/scRNA-Seq/HP4_outs/HP4_data.txt",infile = "txt", outfile = "txt", out_dir="E:/scRNA-Seq/HP4_outs/scImpute_k5/",labeled = FALSE,drop_thre =0.5, Kcluster =5,labels = NULL,ncores = 1)

#cyclone
library(scran)
library(org.Hs.eg.db)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
HSMM_expr_matrix <- read.table("E:/scRNA-Seq/HP4_outs/scImpute_k5/scImpute/scimpute_count.txt")
HSMM_expr_matrix_ok=HSMM_expr_matrix
dim(HSMM_expr_matrix_ok)
genes=read.table("E:/scRNA-Seq/genes.tsv",header=T)
HSMM_expr_matrix_ok$genename=rownames(HSMM_expr_matrix_ok)
HSMM_expr_matrix_ok=merge(genes,HSMM_expr_matrix_ok,by="genename")
dim(HSMM_expr_matrix_ok)
rownames(HSMM_expr_matrix_ok)=HSMM_expr_matrix_ok$ensembl
HSMM_expr_matrix_ok=subset.data.frame(HSMM_expr_matrix_ok,select = -c(ensembl,genename))
dim(HSMM_expr_matrix_ok)
assigned<- cyclone(as.matrix(HSMM_expr_matrix_ok), hs.pairs, gene.names=rownames(HSMM_expr_matrix_ok))
write.table(assigned,"assigned.txt",sep="\t")
assigned=read.delim("assigned.txt",sep="\t",header=T,row.names=1) 
cells=colnames(HSMM_expr_matrix)
rownames(assigned)=cells
write.table(assigned,"cell_sample_sheet.txt",sep="\t") 

col <- character(ncol(HSMM_expr_matrix_ok))
is.G1=grep("G1",assigned$phases)
is.G2M=grep("G2M",assigned$phases)
is.S=grep("S",assigned$phases)
col[is.G1] <- "red"
col[is.G2M] <- "blue"
col[is.S] <- "darkgreen"
plot(assigned$score$G1, assigned$score$G2M, col=col, pch=16)

#seurat
library(Seurat)
setwd("E:/scRNA-Seq/HP4_outs/scImpute_k5")
HSMM_expr_matrix <- read.table("E:/scRNA-Seq/HP4_outs/scImpute_k5/scImpute/scimpute_count.txt")
cell_sample_sheet <- read.table("E:/scRNA-Seq/HP4_outs/scImpute_k5/cell_sample_sheet.txt",header=T,row.names = 1,sep="\t")
HP4<-CreateSeuratObject(raw.data=HSMM_expr_matrix ,min.cells=10,min.genes=200,project="HP4",is.expr=0,meta.data=cell_sample_sheet)
mito.genes<-grep(pattern="^MT-",x=rownames(x=HP4@data),value=TRUE)
percent.mito<-Matrix::colSums(HP4@raw.data[mito.genes,])/Matrix::colSums(HP4@raw.data)
HP4<-AddMetaData(HP4,metadata=percent.mito,col.name="percent.mito")
Total_mRNAs <- Matrix::colSums(HP4@raw.data)
HP4<-AddMetaData(HP4,metadata=Total_mRNAs,col.name="Total_mRNAs")
upper_bound <- 10^(mean(log10(HP4@meta.data$Total_mRNAs)) +2*sd(log10(HP4@meta.data$Total_mRNAs)))
lower_bound <- 10^(mean(log10(HP4@meta.data$Total_mRNAs)) -2*sd(log10(HP4@meta.data$Total_mRNAs)))
HP4<-FilterCells(HP4,subset.names=c("percent.mito","Total_mRNAs"),low.thresholds=c(-Inf,lower_bound),high.thresholds=c(0.05,upper_bound))
HP4<-NormalizeData(HP4,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP4<-AddMetaData(HP4,metadata=HP4@data["HBB",],col.name="HBB_exprs")
HP4<-AddMetaData(HP4,metadata=HP4@data["HBG2",],col.name="HBG2_exprs")
HP4<-AddMetaData(HP4,metadata=HP4@data["HBG1",],col.name="HBG1_exprs")
HP4<-AddMetaData(HP4,metadata=HP4@data["HBA2",],col.name="HBA2_exprs")
HP4<-AddMetaData(HP4,metadata=HP4@data["HBA1",],col.name="HBA1_exprs")
HP4<-AddMetaData(HP4,metadata=HP4@data["GYPA",],col.name="GYPA_exprs")
HP4<-ScaleData(HP4,vars.to.regress="phases")

#HP4<-FindVariableGenes(HP4, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#HP4<-RunPCA(HP4,pc.genes=HP4@var.genes,do.print=T,pcs.print=0,genes.print=0)   
#HP4 =RunTSNE(HP4, do.fast = TRUE,dims.use = 1:10)
#write.table(HP4@dr$tsne@cell.embeddings,"HP4.matrix.txt",sep="\t")
#write.table(as.matrix(HP4@data),"HP4.normalized_expression.txt",sep="\t")
#write.table(as.matrix(HP4@scale.data),"HP4.filtered_expression.txt",sep="\t")

write.table(HP4@meta.data,"HP4.meta.data.txt",sep="\t")
save(HP4,file="HP4.Rdata")

plot(HP4@meta.data$HBG2_exprs,HP4@meta.data$HBB_exprs,pch=20)
plot(HP4@meta.data$HBG1_exprs,HP4@meta.data$HBB_exprs,pch=20)
plot(HP4@meta.data$HBA2_exprs,HP4@meta.data$HBB_exprs,pch=20)
plot(HP4@meta.data$HBA1_exprs,HP4@meta.data$HBB_exprs,pch=20)
plot(HP4@meta.data$GYPA_exprs,HP4@meta.data$HBB_exprs,pch=20)

HP4_N_cells=rownames(HP4@meta.data[which(HP4@meta.data$HBB_exprs <= 6 & HP4@meta.data$HBG2_exprs <= 3 & HP4@meta.data$HBG1_exprs <= 3 & HP4@meta.data$HBA2_exprs <= 5  & HP4@meta.data$HBA1_exprs <= 5 & HP4@meta.data$GYPA_exprs <=3 ),])
HP4_P_cells=setdiff(HP4@cell.names,HP4_N_cells) 
HP4P=SubsetData(HP4,cells.use=HP4_P_cells,subset.raw=T)
HP4P<-NormalizeData(HP4P,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP4P<-ScaleData(HP4P,vars.to.regress="phases")
HP4P<-FindVariableGenes(HP4P, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#HP4P<-RunPCA(HP4P,pc.genes=HP4P@var.genes,do.print=T,pcs.print=0,genes.print=0)   
#HP4P =RunTSNE(HP4P, do.fast = TRUE,dims.use = 1:10)
#write.table(HP4P@dr$tsne@cell.embeddings,"HP4P.matrix.txt",sep="\t")
#write.table(as.matrix(HP4P@data),"HP4P.normalized_expression.txt",sep="\t")
#write.table(as.matrix(HP4P@scale.data),"HP4P.filtered_expression.txt",sep="\t")

write.table(HP4P@meta.data,"HP4P.meta.data.txt",sep="\t")
save(HP4P,file="HP4P.Rdata")

#monocle part
library(monocle)
library(dplyr)
HSMM=importCDS(HP4P, import_all =T) 
HSMM <- estimateSizeFactors(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0)
head(fData(HSMM))
head(pData(HSMM))
summary(pData(HSMM)$nGene)
summary(pData(HSMM)$num_genes_expressed)
summary(fData(HSMM)$num_cells_expressed)

qplot(Total_mRNAs, data = pData(HSMM), geom ="density",xlim=c(0,upper_bound)) + geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)

ggplot(pData(HSMM), aes(num_genes_expressed, Total_mRNAs)) + geom_point()  #We could remove the cells with much higher gene (and UMI) counts as they might be doublets 
x <- pData(HSMM)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
summary(x_1)
df <- data.frame(x = x_1)
ggplot(df, aes(x)) +geom_histogram(bins = 200) +geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

genes_exprs=apply(as.matrix(exprs(HSMM)),1,mean)
genes_exprs=as.matrix(genes_exprs)
high_exprs_genes=rownames(as.matrix(genes_exprs[which(genes_exprs>0.1),]))
HSMM <- setOrderingFilter(HSMM, high_exprs_genes)  #column: use_for_ordering
table(fData(HSMM)$use_for_ordering)
HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F)  # norm_method='log'

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim =50,reduction_method = 'tSNE', verbose = T)  #可加residualModelFormulaStr = "~num_genes_expressed",表示要substrat
HSMM <- clusterCells(HSMM, num_clusters =7)
plot_cell_clusters(HSMM, 1, 2,color = "Cluster",cell_size=1.5)

########Constructing Single Cell Trajectories########
diff_test_res <- differentialGeneTest(HSMM[high_exprs_genes,], fullModelFormulaStr = "~Cluster") 
ordering_genes_diff<- row.names (subset(diff_test_res, qval < 0.1))   #Select genes that are significant at an FDR < 10%
HSMM <- setOrderingFilter(HSMM, ordering_genes_diff)
HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "Cluster")
save(HSMM,file="HSMM.Rdata")

###cellrouter part
write.table(as.matrix(HP4P@data),"HP4P.normalized_expression.txt",sep="\t")
write.table(colnames(HP4P@data),"HP4P.cell_names.csv",sep="\t")
write.table(rownames(HP4P@data),"HP4P.gene_names.csv",sep="\t")

##########some plot start##############
pdf("HP4P.pseudotime_heatmap.pdf")
heatmap=plot_pseudotime_heatmap(HSMM[high_exprs_genes,],num_clusters = 4,cores = 1,show_rownames = T,return_heatmap = T)
dev.off()

c1 <- as.data.frame(cutree(heatmap$tree_row, k=4))
colnames(c1) <- "Cluster"
c1$Gene <- rownames(c1)
write.table(c1,"HP4P.genes_in_heatmap_clusters.txt",sep="\t")

write.table(pData(HSMM),"HP4P.pheno.txt",sep="\t")
write.table(fData(HSMM),"HP4P.feature.txt",sep="\t")

pdf("HP4P.trajactory_by_cluster.pdf")
plot_cell_trajectory(HSMM, color_by = "Cluster")
dev.off()
pdf("HP4P.trajactory_by_Pseudotime.pdf")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()
pdf("HP4P.tSNE.pdf")
plot_cell_clusters(HSMM, 1, 2,color = "Cluster")
dev.off()
pdf("HP4P.6genes.pdf")
plot_cell_clusters(HSMM, 1, 2, color="Cluster", markers=c("SLC4A1","TRFC","EPOR","GYPA","BCAM","ITGA4"))
dev.off()
pdf("HP4P.8genes.pdf")
plot_cell_clusters(HSMM, 1, 2, color="Cluster", markers=c("HK1","NDEL1","EPB41","LGALS9","NEK1","AXIN1","MARK2","USO1","NFE2L1"))
dev.off()
pdf("HP4P.hbgenes.1.pdf")
plot_cell_clusters(HSMM, 1, 2, color="Cluster", markers=c("HBB","HBA2","HBA1","HBG2","HBG1"))
dev.off()
pdf("HP4P.hbgenes.2.pdf")
plot_cell_clusters(HSMM, 1, 2, color="Cluster", markers=c("HBD","HBE1","HBZ","HBM","HBQ1"))
dev.off()
pdf("HP4P.Pseudotime.pdf")
plot_cell_clusters(HSMM, 1, 2,color = "Pseudotime")
dev.off()
pdf("HP4P.cell_cycle.pdf")
plot_cell_clusters(HSMM, 1, 2,color = "phases")
dev.off()

##########plot end##############