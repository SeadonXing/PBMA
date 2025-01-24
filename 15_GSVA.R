# Some required R packages
library(GSVA)
library(msigdbr)
library(GSEABase)
library(limma)

# Step1 : Load Seurat object
scObject <- readRDS("Malig.merge.Seurat.rds")

scObject$TumorGroup_Disease <- paste0(scObject$TumorGroup,"_",scObject$Disease)

# Step2 : AverageExpression
Idents(scObject) <- scObject$TumorGroup_Disease

expr <- AverageExpression(scObject, assays = "RNAmerge", slot = "data")$RNA

expr <- expr[rowSums(expr)>0, ]

expr <- as.matrix(log2(expr+1))

# Step3 : GSVA
H50 <- msigdbr(species = "Homo sapiens", category = "H")

geneset <- split(x = H50$gene_symbol, f = H50$gs_name)

names(geneset) <- sapply(names(geneset), function(x) stringr::str_to_title( gsub("HALLMARK_","",x) ))

GSVAresults <- gsva( expr , geneset , min.sz=5, max.sz=500, method="gsva", kcdf="Gaussian", mx.diff=TRUE, abs.ranking=FALSE, verbose=TRUE, parallel.sz=1)

# Step4 : heatmap
pheatmap::pheatmap( 
	GSVAresults, 
	color = colorRampPalette(c("steelblue", "white", "firebrick"))(100), 
	scale = "row", 
	cluster_rows=TRUE, 
	cluster_cols=FALSE, 
	gaps_col=6,
	clustering_method="ward.D", 
	border_color="black", 
	fontsize_row = 10, 
	fontsize_col = 10, 
	angle_col=90, 
	breaks=seq(-2.5,2.5, length=100)
	)
