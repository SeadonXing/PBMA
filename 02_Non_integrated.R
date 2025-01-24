# Some required R packages
library(Seurat)

# Parallel Computing
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 100 * 1024^3)
options(future.rng.onMisuse = "ignore")

# Step1 : Load seurat object
scObject <- readRDS("BrM.Seurat.rds")
scObject

# Step2 : Perform standard preprocessing (log-normalization), identify variable features)
scObject <- NormalizeData(scObject)

scObject <- FindVariableFeatures(scObject)

# Step3 : Dimensional reduction and clustering
# (1) Using RNA assay
DefaultAssay(object = scObject) <- "RNA"

# (2) Scaling the data and removing unwanted sources of variation
scObject <- ScaleData(object = scObject, model.use="linear")

# (3) Perform linear dimensional reduction
scObject <- RunPCA(object = scObject, npcs = 50)

# (4) Determine statistically significant principal components
ElbowPlot(object = scObject, ndims = 50, reduction = "pca")

# (5) Cluster the cells
scObject <- FindNeighbors(object = scObject, reduction = "pca", dims = 1:30, force.recalc = TRUE)

scObject <- FindClusters(object = scObject, resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))

# (6) Run Non-linear dimensional reduction (UMAP/tSNE)
scObject <- RunTSNE(object = scObject, reduction = "pca", dims = 1:30, check_duplicates = FALSE)

scObject <- RunUMAP(object = scObject, reduction = "pca", dims = 1:30, return.model = TRUE)

# Step4 : Finding differentially expressed genes (cluster biomarkers)
DefaultAssay(object = scObject) <- "RNA"

Idents(scObject) <- scObject$"RNA_snn_res.0.8"

scObject.markers <- FindAllMarkers(object = scObject, 
	only.pos = TRUE, 
	min.pct = 0.2, 
	logfc.threshold = 0.25, 
	max.cells.per.ident=200, 
	assay="RNA", 
	slot = "data")

saveRDS(scObject.markers, file = "markers.rds")
