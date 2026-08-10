# Author: Francesca Petralia
# Affiliation: Icahn School of Medicine at Mount Sinai

library(BayesDeBulk)
library(stringr)

resolve_script_dir <- function(script_name) {
  for (frame in rev(sys.frames())) {
    ofile <- frame$ofile
    if (!is.null(ofile) && basename(ofile) == script_name && file.exists(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)))
  }
  if (file.exists(file.path(getwd(), script_name))) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  stop("Cannot determine script directory for `", script_name, "`. Run from the script folder or with Rscript.", call. = FALSE)
}
script_dir <- resolve_script_dir("Figure7_BayesDeBulk_run.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -- load  Data
load(file.path(script_dir, "Data.rda"))

geneID<-pro[,1]
data<-pro[,-seq(1,4)]

index<-(rowSums(is.na(data))==0)
data<-data[index,]
geneID<-geneID[index]

replicates <- function(data,geneID){ 
  
  # -- When gene/protein replicates exist choose the one with highest interquantile range
  # -- data : (p x n) matrix with rowmanes indicating geneID/proteinID
  
  data<-data[!is.na(geneID),]; geneID<-geneID[!is.na(geneID)];
  dup<-duplicated(geneID)  
  geneDup<-unique(geneID[dup])
  dataNew<-data[is.na(match(geneID,geneDup)),]
  geneNew<-geneID[is.na(match(geneID,geneDup))]
  rownames(dataNew)<-geneNew
  
  
  for (j in 1:length(geneDup)){
    data.dup<-data[geneID==geneDup[j],]
    iqr<-apply(data.dup,1,function(x) {quantile(x,.75)-quantile(x,.25)})
    dataNew<-rbind(dataNew,data.dup[iqr==max(iqr),])
  }
  rownames(dataNew)<-c(geneNew,geneDup)
  
  return(dataNew)
}
pro<-replicates(data,geneID)


gene.rna<-rna[,1]
rna<-rna[,-seq(1,3)]
rownames(rna)<-gene.rna

# -- match RNA and proteome data
mg<-match(colnames(pro),colnames(rna))
rna.new<-rna[,mg[!is.na(mg)]]
pro.new<-pro[,!is.na(mg)]

mg<-match(colnames(pro),colnames(rna))
pro.new<-cbind(pro.new,pro[,is.na(mg)])
rna.new<-cbind(rna.new,matrix(NA,dim(rna.new)[1],sum(is.na(mg))))

colnames(rna.new)<-colnames(pro.new)

mg<-match(colnames(rna),colnames(pro))
rna.new<-cbind(rna.new,rna[,is.na(mg)])
pro.new<-cbind(pro.new,matrix(NA,dim(pro.new)[1],sum(is.na(mg))))
colnames(pro.new)<-colnames(rna.new)

rna<-rna.new
data<-pro.new
sampleid<-colnames(rna)
pro.id<-rownames(data)
gene.id<-rownames(rna)

# -- load cell type signatures for deconvolution

load(file = file.path(script_dir, "Gene_signature.rda"))

cell.type<-unique(c(index.matrix[,1],index.matrix[,2]))

# -- run bayesdebulk
n.iter<-5000
burn.in<-1000
k.fix<-length(cell.type)
data<-t(apply(data,1,function(x)(x-mean(x[!is.na(x)]))/sd(x[!is.na(x)])))
rna<-t(apply(rna,1,function(x)(x-mean(x[!is.na(x)]))/sd(x[!is.na(x)])))

gibbs<-BayesDeBulk(n.iter=n.iter,burn.in=burn.in,Y=list(data,rna),markers=index.matrix,prior=NULL) 

pi.post<-gibbs[[1]]

saveRDS(pi.post, file.path(output_dir, "Figure7_BayesDeBulk_pi_post.rds"))
