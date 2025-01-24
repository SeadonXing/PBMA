# Some required R packages
library(numbat)
library(dplyr)
library(Seurat)
library(ggplot2)
library(glue)
library(data.table)
library(ggtree)
library(stringr)
library(tidygraph)
library(patchwork)

# Step1 : Reload seurat object
scObject <- readRDS("BrM.Seurat.rds")
scObject

# Step2 : subset
DefaultAssay(object = scObject) <- "RNA"

notum <- subset(scObject, FinalCellType %in% c("T cell", "NK", "B cell", "Myeloid"))

# Step3 : Build ref_internal
ref_mat <- notum@assays$RNA@counts
# ref_mat is a gene x cell raw count matrices

cell_annot <- data.frame(cell= colnames(notum),group= as.character(notum$FinalCellType))
# cell_annot is a dataframe with columns "cell" and "group"

ref_internal <- aggregate_counts(ref_mat, cell_annot)

saveRDS(ref_internal, "numbat.ref_internal.rds")
