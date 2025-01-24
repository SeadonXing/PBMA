# Code adapted from https://doi.org/10.1038/s41588-023-01357-3

# Some required R packages
library(Seurat)
library(dplyr)
library(Matrix)
library(NMF)

# Parallel Computing
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 50 * 1024^3)
options(future.rng.onMisuse = "ignore")

# Inputs
pat <- commandArgs(TRUE)[1]
num <- as.numeric(commandArgs(TRUE)[2])

# Step1 : Reload and Pre-filter cells
print( paste0("Start run of : ",pat) )

## malignant
expr_tumor <- readRDS( paste0("BrM.", pat, ".scObject.rds") )

cells <- rownames(metadata)[ metadata$malignant_overlap=="Both_malignant" ]

expr_tumor_pure <- expr_tumor$RNA@counts[ , cells ]

# Step2 : run NMF
nmf_programs <- function(cpm, is.log=F, rank, method="snmf/r", seed=1) {
	if(is.log==F) CP100K_log <- log2(t(t(cpm)/colSums(cpm))*100000+1) else CP100K_log <- cpm
	print(dim(CP100K_log))
	CP100K_log <- CP100K_log[apply(CP100K_log, 1, function(x) length(which(x > 3.5)) > ncol(CP100K_log)*0.02),]
	print(dim(CP100K_log))
	CP100K_log <- CP100K_log - rowMeans(CP100K_log)
	mat.center <- CP100K_log
	CP100K_log[CP100K_log < 0] <- 0
	nmf_programs_obj <- nmf(CP100K_log, rank=rank, method=method, seed=seed)
	return( list(obj=nmf_programs_obj, mat=mat.center) )
}

###
res <- nmf_programs( as.matrix( expr_tumor_pure ), is.log=F, rank=num )
saveRDS( res$obj, paste0(pat,".rank",num,".NMF.rds") )
saveRDS( res$mat, paste0(pat,".rank",num,".Matrix.rds") )

# Step3 : run DEGs
meta.nmf <- function(obj, matrix, k=10, size=10, geneno=100, name="x",...){
  n1<-obj
  
  #Extract top genes per factor
  if(!is.null(geneno)){
    fs <- NMF::extractFeatures(n1, method = geneno)
    names(fs)<-seq(1,k)}else{
      fs<-NMF::extractFeatures(n1)
      names(fs)<-seq(1,k)
      f2<-lapply(fs,length)
      fs<-fs[f2>=size]
    }
  
  #Assign each cell to a factor
  cs<- apply(NMF::coefficients(n1), 2, which.max)
  cs<-cs[cs %in% names(fs)]
  
  #Remove factors with few cells
  cx<-table(cs)>=size
  cx<-names(cx[cx==TRUE])
  cs<-cs[cs %in% cx]
  fs<-fs[cx]
  
  #Assign genes to cells
  jf<-lapply(names(fs),function(x){
    Gene<-rownames(NMF::basis(n1))[fs[[x]]]
    Cluster<-x
    cbind.data.frame(Gene,Cluster)
  })
  jf<-do.call(rbind.data.frame,jf)
  
  #Rename clusters
  cu2<-sort(unique(cs))
  for(i in 1:length(cu2)){
    cs[cs==cu2[i]]<-paste0(name,"_Cluster_",LETTERS[i])
    jf[jf$Cluster==cu2[i],"Cluster"]<-paste0(name,"_Cluster_",LETTERS[i])
  }
  
  #Set logFC
  jf$logFC<-unlist(lapply(unique(jf$Cluster),function(x){
    m2<-matrix[jf[jf$Cluster==x,"Gene"],]
    g1 <- names(cs[cs==x])
    g2 <- names(cs[cs!=x])
    l1<- rowMeans(m2[, g1]) - rowMeans(m2[, g2])
    
  }))

  out<-list(clusters=cs,clustermat=jf)
  return(out)
}

###
ress <- meta.nmf(obj=res$obj, matrix=res$mat, name=pat, k=num, size=5, geneno=100)
saveRDS( ress, paste0(pat,".rank",num,".clustermat.rds") )
print( paste0("Finish for : ", pat, " Rank=", num) )
