# Seurat version 2.3.4
#use scImpute to impute the dropouts before downstream analysis.Assumed that all R packages needed are installed successfully.
library(scImpute)
library("Seurat", lib.loc="D:/software/R-3.5.1/library")
setwd("E:/scRNA-Seq/HP3_outs/scImpute_k5")
HP3.data<-Read10X(data.dir="E:/scRNA-Seq/HP3_outs/outs/filtered_gene_bc_matrices/hg19")
write.table(as.matrix(HP3.data),"E:/scRNA-Seq/HP3_outs/HP3_data.txt",sep="\t")
scimpute("E:/scRNA-Seq/HP3_outs/HP3_data.txt",infile = "txt", outfile = "txt", out_dir="E:/scRNA-Seq/HP3_outs/scImpute_k5/",labeled = FALSE,drop_thre =0.5, Kcluster =5,labels = NULL,ncores = 1)

#cyclone to evaluate the stage of cells
library(scran)
library(org.Hs.eg.db)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
HSMM_expr_matrix <- read.table("E:/scRNA-Seq/HP3_outs/scImpute_k5/scImpute/scimpute_count.txt")
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
write.table(assigned,"cell_sample_sheet.txt",sep="\t",row.names=F)

col <- character(ncol(HSMM_expr_matrix_ok))
is.G1=grep("G1",assigned$phases)
is.G2M=grep("G2M",assigned$phases)
is.S=grep("S",assigned$phases)
col[is.G1] <- "red"
col[is.G2M] <- "blue"
col[is.S] <- "darkgreen"
plot(assigned$score$G1, assigned$score$G2M, col=col, pch=16)

#Data quality control with Seurat package
library(Seurat)
setwd("E:/scRNA-Seq/HP3_outs/scImpute_k5")
HSMM_expr_matrix <- read.table("E:/scRNA-Seq/HP3_outs/scImpute_k5/scImpute/scimpute_count.txt")
cell_sample_sheet <- read.table("E:/scRNA-Seq/HP3_outs/scImpute_k5/cell_sample_sheet.txt",header=T,row.names = 1,sep="\t")
HP3<-CreateSeuratObject(raw.data=HSMM_expr_matrix ,min.cells=10,min.genes=200,project="HP3",is.expr=0,meta.data=cell_sample_sheet)
mito.genes<-grep(pattern="^MT-",x=rownames(x=HP3@data),value=TRUE)
percent.mito<-Matrix::colSums(HP3@raw.data[mito.genes,])/Matrix::colSums(HP3@raw.data)
HP3<-AddMetaData(HP3,metadata=percent.mito,col.name="percent.mito")
Total_mRNAs <- Matrix::colSums(HP3@raw.data)
HP3<-AddMetaData(HP3,metadata=Total_mRNAs,col.name="Total_mRNAs")
upper_bound <- 10^(mean(log10(HP3@meta.data$Total_mRNAs)) +2*sd(log10(HP3@meta.data$Total_mRNAs)))
lower_bound <- 10^(mean(log10(HP3@meta.data$Total_mRNAs)) -2*sd(log10(HP3@meta.data$Total_mRNAs)))
HP3<-FilterCells(HP3,subset.names=c("percent.mito","Total_mRNAs"),low.thresholds=c(-Inf,lower_bound),high.thresholds=c(0.05,upper_bound))
HP3<-NormalizeData(HP3,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP3<-AddMetaData(HP3,metadata=HP3@data["HBB",],col.name="HBB_exprs")
HP3<-AddMetaData(HP3,metadata=HP3@data["HBG2",],col.name="HBG2_exprs")
HP3<-AddMetaData(HP3,metadata=HP3@data["HBG1",],col.name="HBG1_exprs")
HP3<-AddMetaData(HP3,metadata=HP3@data["HBA2",],col.name="HBA2_exprs")
HP3<-AddMetaData(HP3,metadata=HP3@data["HBA1",],col.name="HBA1_exprs")
HP3<-AddMetaData(HP3,metadata=HP3@data["GYPA",],col.name="GYPA_exprs")
HP3<-ScaleData(HP3,vars.to.regress="phases")

write.table(HP3@meta.data,"HP3.meta.data.txt",sep="\t")
save(HP3,file="HP3.Rdata")

plot(HP3@meta.data$HBG2_exprs,HP3@meta.data$HBB_exprs,pch=20)
plot(HP3@meta.data$HBG1_exprs,HP3@meta.data$HBB_exprs,pch=20)
plot(HP3@meta.data$HBA2_exprs,HP3@meta.data$HBB_exprs,pch=20)
plot(HP3@meta.data$HBA1_exprs,HP3@meta.data$HBB_exprs,pch=20)
plot(HP3@meta.data$GYPA_exprs,HP3@meta.data$HBB_exprs,pch=20)

#N:globin-Negtive cells.P:globin-Positive cells
HP3_N_cells=rownames(HP3@meta.data[which(HP3@meta.data$HBB_exprs <= 5 & HP3@meta.data$HBG2_exprs <= 6 & HP3@meta.data$HBG1_exprs <= 6 & HP3@meta.data$HBA2_exprs <= 5  & HP3@meta.data$HBA1_exprs <= 5 & HP3@meta.data$GYPA_exprs ==0 ),])
HP3_P_cells=setdiff(HP3@cell.names,HP3_N_cells)
HP3P=SubsetData(HP3,cells.use=HP3_P_cells,subset.raw=T)
HP3P<-NormalizeData(HP3P,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP3P<-ScaleData(HP3P,vars.to.regress="phases")
HP3P<-FindVariableGenes(HP3P, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

write.table(HP3P@meta.data,"HP3P.meta.data.txt",sep="\t")
save(HP3P,file="HP3P.Rdata")

###prepare normalized expression for Cellrouter input
write.table(as.matrix(HP3P@data),"HP3P.normalized_expression.txt",sep="\t")
write.table(colnames(HP3P@data),"HP3P.cell_names.csv",sep="\t")
write.table(rownames(HP3P@data),"HP3P.gene_names.csv",sep="\t")
