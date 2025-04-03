# Some required R packages
library(Seurat)
library(philentropy)
library(ggplot2)

# some functions
'%notin%' <- Negate('%in%')
get_field <- function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])

# Step1 : Calculate M1/M2 state scores
scObject <- scObject <- readRDS("Myeloid.merge.Seurat.rds")

M1 <- list( c("IL23A", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86", "IL1A", "IL1B", "IL6", "CCL5", "IRF5", "IRF1", "CD40", "IDO1", "KYNU", "CCR7") )
scObject <- AddModuleScore(object = scObject, features = M1, name = 'M1_score')

M2 <- list( c("IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", "CCL24", "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", "CTSA", "CTSB", "CTSC", "CTSD", 
	"TGFB1", "TGFB2", "TGFB3", "MMP14", "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASLG", "TNFSF12", "TNFSF8", "CD276", "VTCN1", "MSR1", "FN1", "IRF4", "MRC1", "CD163") )
scObject <- AddModuleScore(object = scObject, features = M2, name = 'M2_score')

# Step2 : Capture differences in cellular states between BrM and PTs across various cancer types
whichstat <- "M1_score1"
clusters <- c("Mono-Resting","Mono-Inflammation","Mono-IFN","TAMacro-Hypoxic","TAMacro-Lipid","TAMacro-Metal","TAMacro-IFNI",
			  "TAMacro-IFNII","TAMacro-Inflammation","TAMacro-Resident","TAMacro-Phagocytosis","Cycling-G1S","Cycling-G2M")

vln_df <- scObject@meta.data %>% mutate(gene=scObject@meta.data[,whichstat]) 
	%>% filter(FinalCellType%in%clusters, TumorGroup%notin%c("Others","Sarcoma")) 
	%>% rename("Disease"="Feature")
vln_df$Disease_Group_ID <- paste0(vln_df$Feature, "#", vln_df$TumorGroup, "#", vln_df$SampleID)

plotdat <- split(vln_df$gene,vln_df$Disease_Group_ID) %>% lapply(.,function(x) quantile(as.numeric(x), seq(0.05,0.95,0.05)))
	%>% do.call(rbind,.) %>% melt() %>% as.data.frame()
plotdat$disease <- sapply( as.character(plotdat$Var1), get_field, 1, "#")
plotdat$group <- sapply( as.character(plotdat$Var1), get_field, 2, "#")
plotdat$group <- factor(plotdat$group, levels=c("LC","BC","MEL","CRC","ESCC","HCC"))

p1 <- ggplot(plotdat, aes(x = disease, y = value)) + 
	 geom_violin(aes(fill=disease), trim=TRUE, scale = "width", adjust=1) + 
	 geom_boxplot(aes(fill=disease), outlier.shape=NA, color="black", width=0.3) +
	 theme_classic() + scale_fill_manual() +
	 facet_wrap( ~ group, nrow=1, scale="fixed") +
	 stat_compare_means(method="wilcox.test", label = "p.format")+
	 theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1, size=10))

# Step3 : Calculate Transformation Index (TI) scores
metadat <- scObject@meta.data %>% dplyr::select(Disease,TumorGroup,FinalCellType,M1_score1,M2_score1) %>%
	mutate(M1_scalescore=scales::rescale(M1_score1, to = c(0, 1)), M2_scalescore=scales::rescale(M2_score1, to = c(0, 1))) # rescale to a range of 0-1

## data
vardata <- metadat %>% group_by(Disease,TumorGroup,FinalCellType) 
	%>% summarise(M1inflam=mean(M1_scalescore), M2anti=mean(M2_scalescore)) %>% split(.,.[,"Disease"])

whichstat <- "M1inflam"
dat01.stat <- vardata$BrM %>% reshape2::dcast( FinalCellType ~ TumorGroup, value.var=whichstat) %>% dplyr::select(-Sarcoma,-Others) %>% column_to_rownames("FinalCellType")
colnames(dat01.stat) <- paste0("BrM_",colnames(dat01.stat))
dat02.stat <- vardata$TIS %>% reshape2::dcast( FinalCellType ~ TumorGroup, value.var=whichstat) %>% column_to_rownames("FinalCellType")
colnames(dat02.stat) <- paste0("PT_",colnames(dat02.stat))
identical(rownames(dat01.stat),rownames(dat02.stat)) # Should be TRUE !!!
dat <- as.data.frame( t(cbind(dat01.stat, dat02.stat)) ) %>% dplyr::select(all_of(clusters))
dat
is.na(dat) # check NA value

## sign
signs <- apply(dat,1,median)[1:6]-apply(dat,1,median)[7:12]
signs <- ifelse(signs>=0,1,-1)

## distance
distval <- philentropy::distance(dat, method = "soergel")
rownames(distval) = colnames(distval) <- rownames(dat)

## plot data
bardat <- melt(distval)[c(7,20,33,46,59,72),]
bardat$TumorGroup <- sapply(as.character(bardat$Var2), get_field, 2, "_")
bardat$signs <- signs
bardat$TI_score <- (bardat$signs)*(bardat$value)
p2 <- ggplot(bardat, aes(x=factor(Var2,levels=bardat$Var2[order(bardat$TI_score)]),y=TI_score)) + 
  geom_bar(aes(fill = TumorGroup),stat="identity") + scale_fill_manual() + 
  theme(axis.line = element_line(colour = "black"), axis.text.x=element_text(angle = 90)) + coord_flip() +
  theme_bw()+labs(x="", y="TI score", fill="")

## final plots
p2+p1
