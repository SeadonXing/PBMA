# Some required R packages
library(Seurat)
library(SeuratObject)
library(harmony)

# Parallel Computing
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 50 * 1024^3)
options(future.rng.onMisuse = "ignore")

# Step1 : Load seurat object
scObject.integrated <- readRDS("BrM.Seurat.rds")
scObject.integrated

# Step2 : Perform standard preprocessing (log-normalization, identify variable features, scale and PCA)
scObject.integrated <- NormalizeData(scObject.integrated) %>% 
	FindVariableFeatures() %>% 
	ScaleData() %>% 
	RunPCA(npcs = 30, verbose = FALSE)

# Step3 : Dimensional reduction and clustering
set.seed(123)

scObject.integrated <- RunHarmony(scObject.integrated, 
	group.by.vars = "SampleID", 
	reduction = "pca", 
	dims.use = 1:30, 
	assay.use = "RNA", 
	project.dim = TRUE,
	max.iter.harmony = 10, 
	epsilon.cluster=-Inf, 
	epsilon.harmony=-Inf, 
	plot_convergence=TRUE)

scObject.integrated <- RunUMAP(scObject.integrated, reduction = "harmony", dims = 1:20, return.model = TRUE)

scObject.integrated <- FindNeighbors(scObject.integrated, reduction = "harmony", dims = 1:20)

scObject.integrated <- FindClusters(scObject.integrated, resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))

# Step4 : Finding differentially expressed genes (cluster biomarkers)
DefaultAssay(object = scObject.integrated) <- "RNA"

Idents(scObject.integrated) <- scObject.integrated$"RNA_snn_res.0.5"

scObject.integrated.markers <- FindAllMarkers(object = scObject.integrated, 
	only.pos = TRUE, 
	min.pct = 0.2, 
	logfc.threshold = 0.25, 
	max.cells.per.ident=200, 
	assay="RNA", 
	slot = "data")

saveRDS(scObject.integrated.markers, file = "markers.rds")
