# Some required R packages
library(Seurat)
library(dplyr)
library(stringr)
library(gplots)
library(ggplot2)
library(scales)
library(viridis)
library(DoubletFinder)
library(patchwork)

pat <- commandArgs(TRUE)[1]

# Step1 : PrepareFiles
scObject <- readRDS( paste0("BrM.",pat,".scObject.rds") )

seurat <- CreateSeuratObject(counts = scObject$RNA@counts, project = pat)

seurat <- AddMetaData(seurat, scObject$FinalCellType, col.name="FinalCellType")

# Step2 : Identify and remove doublets
dir.create(pat)

## 1. Pre-process Seurat object (standard)
DefaultAssay(seurat) <- "RNA"

seurat <- NormalizeData(seurat, verbose=FALSE)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose=FALSE)

seurat <- ScaleData(seurat, verbose=FALSE)

seurat <- RunPCA(seurat, verbose=FALSE)

seurat <- RunUMAP(seurat, dims = 1:30, verbose=FALSE)

## 2. pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(seurat, PCs = 1:30, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

bcmvn <- find.pK(sweep.stats)

mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

## 3. Homotypic Doublet Proportion Estimate
DoubletRate <- ncol(seurat)*8*1e-6 ## ~1000 cells, and a multiplet rate of ~0.8% according to 10X Genomics

nExp_poi <- round(DoubletRate*ncol(seurat))

annotations <- seurat$FinalCellType

homotypic.prop <- modelHomotypic(annotations)

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 4. Run DoubletFinder with varying classification stringencies
seurat <- doubletFinder_v3(seurat, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seurat <- doubletFinder_v3(seurat, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

seurat$doubletfinder <- seurat@meta.data[,paste0("DF.classifications_0.25_",mpK,"_",nExp_poi)]

seurat$doubletfinder.adj <- seurat@meta.data[,paste0("DF.classifications_0.25_",mpK,"_",nExp_poi.adj)]

table(seurat$doubletfinder)

table(seurat$doubletfinder.adj)

## 5. Save results
write.table( seurat@meta.data[,"doubletfinder",drop=FALSE], paste0(pat,"/",pat,"_doubletfinder.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t" )

write.table( seurat@meta.data[,"doubletfinder.adj",drop=FALSE], paste0(pat,"/",pat,"_doubletfinder.adj.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t" )
