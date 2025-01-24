# Some required R packages
library(Seurat)

# Parallel Computing
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 100 * 1024^3)
options(future.rng.onMisuse = "ignore")

# Step1 : Integration to remove batch effects : RPCA
scObject <- readRDS("BrM.Seurat.rds")
scObject

# (1) Split Seurat Object by SampleID
scObject.list <- SplitObject(object = scObject, split.by = "SampleID")

# (2) RPCA & Reference
# 1. normalize and identify variable features for each dataset independently
scObject.list <- lapply(X = scObject.list, FUN = function(x) {
	x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})

# 2. select features for integration and run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = scObject.list)
scObject.list <- lapply(X = scObject.list, FUN = function(x) {
	x <- ScaleData(x, features = features, verbose = FALSE)
	x <- RunPCA(x, features = features, verbose = FALSE)
})

# 3. Identify anchors
names(scObject.list)

scObject.anchors <- FindIntegrationAnchors(object.list = scObject.list, anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5, 
	reference = c(95,81,40,   83,2,8,   86,31,88,7,103,   91,58))

# 4. Integration of many cell datasets
scObject.integrated <- IntegrateData(anchorset = scObject.anchors, dims = 1:30)

# Step2 : Dimensional reduction and clustering
# (1) Using integrated assay
DefaultAssay(object = scObject.integrated) <- "integrated"

# (2) Scaling the data and removing unwanted sources of variation
scObject.integrated <- ScaleData(object = scObject.integrated, model.use="linear")

# (3) Perform linear dimensional reduction
scObject.integrated <- RunPCA(object = scObject.integrated, npcs = 50)

# (4) Determine statistically significant principal components
ElbowPlot(object = scObject.integrated, ndims = 50, reduction = "pca")

# (5) Cluster the cells
scObject.integrated <- FindNeighbors(object = scObject.integrated, reduction = "pca", dims = 1:30, force.recalc = TRUE)

scObject.integrated <- FindClusters(object = scObject.integrated, resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))

# (6) Run Non-linear dimensional reduction (UMAP/tSNE)
scObject.integrated <- RunTSNE(object = scObject.integrated, reduction = "pca", dims = 1:30, check_duplicates = FALSE)

scObject.integrated <- RunUMAP(object = scObject.integrated, reduction = "pca", dims = 1:30, return.model = TRUE)

# Step3 : Finding differentially expressed genes (cluster biomarkers)
DefaultAssay(object = scObject.integrated) <- "RNA"

Idents(scObject.integrated) <- scObject.integrated$"integrated_snn_res.0.8"

scObject.integrated.markers <- FindAllMarkers(object = scObject.integrated, 
	only.pos = TRUE, 
	min.pct = 0.2, 
	logfc.threshold = 0.25, 
	max.cells.per.ident=200, 
	assay="RNA", 
	slot = "data")

saveRDS(scObject.integrated.markers, file = "markers.rds")
