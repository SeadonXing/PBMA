# Some required R packages
library(infercnv)
library(Seurat)

pat <- commandArgs(TRUE)[1]

# Step1 : PrepareFiles
scObject <- readRDS( paste0("BrM.",pat,".scObject.rds") )
scObject

scObject$anno <- as.character(scObject$FinalCellType)

scObject$anno[ scObject$anno %in% c("MTCs","Fibroblast","Endothelial","Glial cell") ] <- "Tumor Candidate"

write.table( scObject@meta.data[,"anno",drop=FALSE], paste0(pat,"_annotations_file.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )

table(scObject$FinalCellType)

table(scObject$anno)

ref_group_names <- setdiff(names(table(as.character(scObject$anno))),"Tumor Candidate")

ref_group_names

# Step2 : CreateInfercnvObject
infercnv_obj = CreateInfercnvObject(

	raw_counts_matrix= GetAssayData(scObject, slot="counts"),
	annotations_file= paste0(pat,"_annotations_file.txt"),
	delim="\t",
	gene_order_file= system.file("extdata", "hg38_gencode_v27.txt", package = "infercnv"),
	ref_group_names= ref_group_names,
	max_cells_per_group = NULL,
	min_max_counts_per_cell = c(100, +Inf),
	chr_exclude = c("chrX", "chrY", "chrM")
	)

# Step3 : infercnv
infercnv_obj = infercnv::run(

	infercnv_obj,
	cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10X Genomics
	out_dir= paste0(pat,"_inferCNV"), 
	cluster_by_groups=FALSE, # If observations are defined according to groups (ie. patients), each group of cells will be clustered separately.
	denoise=TRUE,
	HMM=TRUE, # when set to True, runs HMM to predict CNV level (default: FALSE)
	analysis_mode="subclusters", # options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions. default: samples (fastest, but subclusters is ideal)
	tumor_subcluster_partition_method="leiden",
	leiden_resolution = 0.01,
	output_format="pdf",
	num_threads=16
	)
