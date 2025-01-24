# Some required R packages
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(cowplot)
library(ggplot2)

# Parallel Computing
library(future)
availableWorkers()
availableCores()
options(future.globals.maxSize= 30 * 1024^3)
plan(strategy = "multicore", workers = 30)

# Step1 : Define gene set
CIN70 <- c(
    "TPX2", "PRC1", "FOXM1", "CDK1", "RAB5IF", "TGIF2", "MCM2", "H2AFZ", "TOP2A",
    "PCNA", "UBE2C", "MELK", "TRIP13", "NCAPD2", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A",
    "CDC45", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "CKAP5",
    "NUP205", "CDC20", "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "EZH2", "CTPS1",
    "DKC1", "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "NCAPH", "MCM10", "LSM4",
    "NCAPG2", "ASF1B", "ZWINT", "PBK", "ZWILCH", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2",
    "RAD21", "ACTL6A", "GPI", "PDCD2L", "SRSF2", "HDGF", "NXT1", "NEK2", "DHCR7", "AURKA",
    "NDUFAB1", "NEMP1", "KIF4A"
  )

# Step2 : AddModuleScore
geneset <- list("CIN70"=CIN70)

scObject <- AddModuleScore(object = scObject, features = geneset, name = names(geneset))

# Step3 : Plots
v1 <- VlnPlot(object = scObject, features = "CIN70", group.by="Disease", pt.size=0)+
	geom_boxplot(outlier.size=0, color="black", width=0.3)

v2 <- VlnPlot(object = scObject, features = "CIN70", group.by="TumorGroup_Disease", pt.size=0)+
	geom_boxplot(outlier.size=0, color="black", width=0.3)
