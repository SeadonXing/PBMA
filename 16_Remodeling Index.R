# Some required R packages
library(Seurat)
library(philentropy)

# Parallel Computing
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 100 * 1024^3)
options(future.rng.onMisuse = "ignore")

# Step1 : Load metadata
brmmeta <- readRDS("BrM.metadata.rds")
brmmeta$Disease <- "BrM"
ptmeta <- readRDS("PT.metadata.rds")
ptmeta$Disease <- "PT"
cols <- c("Disease","CancerGroup","FinalCellType")
metadat <- rbind(brmmeta[,cols], ptmeta[,cols])

## merge
metadat$RItype <- "Other"
metadat$RItype[ metadat$FinalCellType%in%c("MTC","PTC") ] <- "Tumor"
metadat$RItype[ metadat$FinalCellType%in%c("T_NK","B_cell","Myeloid") ] <- "Immune"
metadat$RItype[ metadat$FinalCellType%in%c("Endothelial","Fibroblast") ] <- "Stromal"

# Step2 : Proportion
dat01 <- metadat %>% filter(RItype!="Other", Disease=="BrM")
dat01 <- as.data.frame.matrix(table(dat01$RItype, dat01$CancerGroup)) %>% dplyr::select(-Sarcoma, -Other)
dat01.prop <- t(t(dat01)/colSums(dat01))*100
colnames(dat01.prop) <- paste0("BrM_",colnames(dat01.prop))

dat02 <- metadat %>% filter(RItype!="Other", Disease=="PT")
dat02 <- as.data.frame.matrix(table(dat02$RItype, dat02$CancerGroup)) %>% dplyr::select(-Sarcoma, -Other)
dat02.prop <- t(t(dat02)/colSums(dat02))*100
colnames(dat02.prop) <- paste0("PT_",colnames(dat02.prop))

# Step3 : Calculate RI
identical(rownames(dat01),rownames(dat02)) # TRUE
dat <- as.data.frame( t(cbind(dat01.prop, dat02.prop)) )

distval <- philentropy::distance(dat, method = "soergel")
rownames(distval) = colnames(distval) <- rownames(dat)
distval

# Step4 : Chi-square test
for(i in c("LC","BC","MEL","CRC","ESCC","HCC")){
	dat <- metadat %>% filter(FinalCellType!="Glial_cell", CancerGroup==i) %>% 
		select(RItype,Disease) %>% table() %>% as.data.frame() %>% group_by(Disease) %>% 
		mutate(Percent = prop.table(Freq) * 100)

	data_matrix <- reshape2::dcast(dat, Disease ~ RItype, value.var="Percent") %>% 
		tibble::column_to_rownames("Disease") %>% round()

	print(chisq.test(data_matrix,simulate.p.value=FALSE))

}
