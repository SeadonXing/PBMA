# Some required R packages
library(Seurat)

# Step1 : Reload scObject
scObject <- readRDS("BrM.Seurat.rds")
scObject

# Step2 : Split Seurat object
table(scObject$SampleID)

DefaultAssay(object = scObject) <- "RNA"

scObject.list <- SplitObject(object = scObject, split.by = "SampleID")

names(scObject.list)

# Step3 : saveRDS
for(x in names(scObject.list)){

	saveRDS(scObject.list[[x]], paste0("BrM.",x,".scObject.rds"))

}
