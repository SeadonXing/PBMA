# Some required R packages
library(Seurat)
library(SeuratWrappers)

# Parallel Computing
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 50 * 1024^3)
options(future.rng.onMisuse = "ignore")

# Step1 : Load RDS
scObject <- readRDS("BrM.Seurat.rds")
scObject

# Step2 : Perform standard preprocessing (log-normalization, identify variable features)
scObject <- NormalizeData(scObject)

scObject <- FindVariableFeatures(scObject)

scObject.list <- SplitObject(scObject, split.by = "SampleID")

names(scObject.list)

# Step3 : Dimensional reduction and clustering
set.seed(123)

scObject.integrated <- RunFastMNN(object.list = scObject.list)

scObject.integrated <- RunUMAP(scObject.integrated, reduction = "mnn", dims = 1:30)

scObject.integrated <- FindNeighbors(scObject.integrated, reduction = "mnn", dims = 1:30)

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
