################################################################################
# Causal kinase analysis
#
# Author: Shrabanti Chowdhury
# Affiliation: Icahn School of Medicine at Mount Sinai
# Purpose: Reproduce Figure 5 and Figure S5 causal kinase analysis panels.
# Output: Panel-labeled PDF files in the Figure 5 analysis folder.
################################################################################

####### libraries to load #######
library(ggplot2)
library(igraph)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(metap)
library(ggrepel)
library(ggpubr)
library(survival)
library(survminer)

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
script_dir <- resolve_script_dir("Figure5_S5_causal_kinase_analysis.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

save_heatmap_pdf <- function(ht, file, width = 6, height = 6) {
  pdf(file.path(output_dir, file), width = width, height = height, useDingbats = TRUE)
  draw(ht)
  dev.off()
}

save_ggplot_pdf <- function(plot, file, width = 6, height = 5) {
  pdf(file.path(output_dir, file), width = width, height = height, useDingbats = TRUE)
  print(plot)
  dev.off()
}

draw_kinase_correlation_heatmap <- function(cor_mat, fdr_mat, file) {
  sig_legend <- Legend(
    labels = "Significant assoc\n(FDR < 0.10)",
    type = "points",
    pch = 16,
    size = unit(3, "mm"),
    legend_gp = gpar(col = "black", fill = "black")
  )

  ht <- Heatmap(
    as.matrix(cor_mat),
    name = "Correlation",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    na_col = "grey75",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 18),
    column_names_gp = gpar(fontsize = 18),
    column_names_rot = 90,
    rect_gp = gpar(col = "white", lwd = 2),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(cor_mat[i, j]) && fdr_mat[i, j] < 0.1) {
        grid.points(
          x, y,
          pch = 16,
          size = unit(2.4, "mm"),
          gp = gpar(fill = "black", col = "black")
        )
      }
    },
    height = nrow(cor_mat) * unit(8, "mm"),
    width = ncol(cor_mat) * unit(8, "mm"),
    column_order = colnames(cor_mat),
    row_order = rownames(cor_mat),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = "Correlation",
      title_gp = gpar(fontsize = 18),
      labels_gp = gpar(fontsize = 14),
      at = c(1, 0.5, 0, -0.5, -1),
      direction = "horizontal",
      legend_width = unit(45, "mm")
    )
  )

  pdf(file.path(output_dir, file), width = 4.3, height = 6.2, useDingbats = TRUE)
  draw(
    ht,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    annotation_legend_list = list(sig_legend),
    padding = unit(c(4, 4, 4, 4), "mm")
  )
  dev.off()
}

####### Figure 5A ##########
load("data/cov_kinase.RData")
ind<- union(which(var.tab$days > (365*1)), intersect(which(var.tab$days < (365*0.5)), which(var.tab$status==1)))

var.tab2<- var.tab[ind,]

age.PED<- ifelse(var.tab2$age.group=="[0,15]", 1, 0)
age.ADOL<- ifelse(var.tab2$age.group=="(15,26]", 1, 0)
age.YA<- ifelse(var.tab2$age.group=="(26,40]", 1, 0)
age.ADULT<- ifelse(var.tab2$age.group=="(40,62]", 1, 0)

var.tab2$SR<- ifelse(var.tab2$days > (365*1), 0, 1)

cortical<- ifelse(var.tab2$location==1,1,0)
midline<- ifelse(var.tab2$location==2,1,0)
other<- ifelse(var.tab2$location==3,1,0)

var.tab2$grade2<- ifelse(var.tab2$grade==3,1,0)

load("data/selected_kinase.RData")
top.kin <- kin.f

load("data/data_kea3.RData")

data2<- data2.sub[,which(colnames(data2.sub) %in% var.tab2$sample_id)]
data3<- data2[which(rownames(data2) %in% kin.f),order(match(colnames(data2), var.tab2$sample_id))]
identical(var.tab2$sample_id, colnames(data3))
# "Sex"=var.tab2$gender, 
datay3<- rbind.data.frame("PED"=age.PED, "ADO"=age.ADOL, "YA"=age.YA, "Sex"=var.tab2$gender, "grade"=var.tab2$grade2, "cortical"=cortical, "midline"=midline, "ATRX"=var.tab2$ATRX, "TP53"=var.tab2$TP53, "H3"=var.tab2$H3F3A, "IDH1"=var.tab2$IDH1, data3, "Prognosis"=var.tab2$SR)

black<- matrix(F,dim(datay3)[1],dim(datay3)[1])
black[8:11, 1:7]<- T # age and sex grade loc <-- (no) mut

black[12:197, c(1:11)]<- T # age, sex, mut, loc <-- (no) kin
black[198, 1:11]<- T # age, sex, mut, kin <-- (no) surv

black[1:7, 1:7]<- T # age and sex and loc <--> (no) age and sex and loc
black[8:11, 8:11]<- T # mut <--> (no) mut

rownames(black)<- colnames(black)<- rownames(datay3) 

if (requireNamespace("dagbagM", quietly = TRUE)) {
  res.hc<- dagbagM::hc_boot_parallel(Y=t(as.matrix(datay3)), n.boot=100, nodeType=c(rep("b",11),rep("c",dim(data3)[1]), rep("b",1)),  whiteList=NULL, blackList=black, standardize=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)

  adj=dagbagM::score_shd(res.hc,threshold=0.4, whitelist = NULL, blacklist = black)$adj.matrix
  rownames(adj)<- colnames(adj)<- rownames(datay3)
} else {
  message("Package `dagbagM` is not installed; using included consensus DAG inputs for Figure 5 network rendering.")
}

load("data/consensus_DAG.RData")

library(igraph)
obj<- graph_from_adjacency_matrix(as.matrix(dag), mode = "directed", weighted = TRUE)

pdf(file.path(output_dir, "Figure5A_causal_kinase_network.pdf"),width = 20, height = 20, useDingbats=T)
plot.igraph(obj,vertex.size=10, edge.length=10, edge.color="black", vertex.label.cex=1.5,vertex.color=V(obj)$color,layout=layout.fruchterman.reingold(obj, niter=10000))
dev.off()

####### Figure 5B, S5B #########
load("data/consensus_DAG.RData")

mod1.kinase<- c("CDK8", "ATM", "ATR", "NTRK2", "GUCY2D", "PLK2", "LCK", "Prognosis")
mod1<- dag[rownames(dag) %in% mod1.kinase, colnames(dag) %in% mod1.kinase]
obj<- graph_from_adjacency_matrix(as.matrix(mod1), mode = "directed", weighted = TRUE)

# %%%%%% 5B %%%%%%%%
pdf(file.path(output_dir, "Figure5B_module_I_causal_kinase_network.pdf"),width = 20, height = 20, useDingbats=T)
plot.igraph(obj,vertex.size=10, edge.length=10, edge.color="black", vertex.label.cex=1.5,vertex.color=V(obj)$color,layout=layout.fruchterman.reingold(obj, niter=10000))
dev.off()
# %%%%%%%%%%%%%%%%%%

mod2.kinase<- c("MAP2K1", "MAP2K2", "JAK2", "RAF1", "FLT1", "PLK2", "GUCY2D", "Prognosis")
mod2<- dag[rownames(dag) %in% mod2.kinase, colnames(dag) %in% mod2.kinase]
obj<- graph_from_adjacency_matrix(as.matrix(mod2), mode = "directed", weighted = TRUE)

# %%%%%%%% S5B %%%%%%%%%%
pdf(file.path(output_dir, "FigureS5B_module_II_causal_kinase_network.pdf"),width = 20, height = 20, useDingbats=T)
plot.igraph(obj,vertex.size=10, edge.length=10, edge.color="black", vertex.label.cex=1.5,vertex.color=V(obj)$color,layout=layout.fruchterman.reingold(obj, niter=10000))
dev.off()
# %%%%%%%%%%%%%%%%%%
####### Figure S5C, D ########
load("data/consensus_DAG_male.RData")

mod1.kinase<- c("CDK8", "ATM", "ATR", "NTRK2", "NTRK3", "GUCY2D", "PLK2", "LCK", "Prognosis")
mod1<- dag[rownames(dag) %in% mod1.kinase, colnames(dag) %in% mod1.kinase]
obj<- graph_from_adjacency_matrix(as.matrix(mod1), mode = "directed", weighted = TRUE)

# %%%%%% S5C %%%%%%%%
pdf(file.path(output_dir, "FigureS5C_male_module_I_causal_kinase_network.pdf"),width = 20, height = 20, useDingbats=T)
plot.igraph(obj,vertex.size=10, edge.length=10, edge.color="black", vertex.label.cex=1.5,vertex.color=V(obj)$color,layout=layout.fruchterman.reingold(obj, niter=10000))
dev.off()
# %%%%%%%%%%%%%%%%%%

mod2.kinase<- c("MAP2K1", "MAP2K2", "JAK2", "RAF1", "FLT1", "PLK2", "GUCY2D", "Prognosis")
mod2<- dag[rownames(dag) %in% mod2.kinase, colnames(dag) %in% mod2.kinase]
obj<- graph_from_adjacency_matrix(as.matrix(mod2), mode = "directed", weighted = TRUE)

# %%%%%%%% S5D %%%%%%%%%%
pdf(file.path(output_dir, "FigureS5D_male_module_II_causal_kinase_network.pdf"),width = 20, height = 20, useDingbats=T)
plot.igraph(obj,vertex.size=10, edge.length=10, edge.color="black", vertex.label.cex=1.5,vertex.color=V(obj)$color,layout=layout.fruchterman.reingold(obj, niter=10000))
dev.off()
# %%%%%%%%%%%%%%%%%%

######## Figure 5C ###########
load("selected_kinase.RData")
xia<- read.csv("data/xiaoyu_corr.csv")
sig.corr<- xia[xia$eFDR_cancer < 0.1, c(1:2,5)]

mod1<- c("ATM", "ATR", "CDK8","PLK2", "GUCY2D", "LCK", "NTRK2")

xia.kin<- xia[intersect(which(xia$kinase %in% top.kin), which(xia$kinase2 %in% top.kin)),]

xia.kin1<- xia.kin[intersect(which(xia.kin$kinase %in% mod1), which(xia.kin$kinase2 %in% mod1)),]

mat1<- xia.kin1
cor1<- matrix(NA, length(mod1), length(mod1))
for(i in 1:length(mod1)){
  for(j in 1:length(mod1)){
    tempi<- which(mat1$kinase == mod1[i])
    tempj<- which(mat1$kinase2 == mod1[j])
    if(length(intersect(tempi, tempj)) != 0){
      cor1[i,j]<- mat1$cor_cancer[intersect(tempi, tempj)]
    }
  }
} 

fdr1<- matrix(NA, length(mod1), length(mod1))
for(i in 1:length(mod1)){
  for(j in 1:length(mod1)){
    tempi<- which(mat1$kinase == mod1[i])
    tempj<- which(mat1$kinase2 == mod1[j])
    if(length(intersect(tempi, tempj)) != 0){
      fdr1[i,j]<- mat1$eFDR_cancer[intersect(tempi, tempj)]
    }
  }
} 

rownames(cor1)<- colnames(cor1)<- rownames(fdr1)<- colnames(fdr1)<- mod1

draw_kinase_correlation_heatmap(cor1, fdr1, "Figure5C_module_I_kinase_correlation_heatmap.pdf")

######## Figure S5E ###########
xia<- read.csv("data/xiaoyu_corr.csv")
sig.corr<- xia[xia$eFDR_cancer < 0.1, c(1:2,5)]

mod1<- c("FLT1", "GUCY2D", "JAK2", "MAP2K1", "MAP2K2", "PLK2",  "RAF1")

xia.kin<- xia[intersect(which(xia$kinase %in% top.kin), which(xia$kinase2 %in% top.kin)),]

xia.kin1<- xia.kin[intersect(which(xia.kin$kinase %in% mod1), which(xia.kin$kinase2 %in% mod1)),]

mat1<- xia.kin1
cor1<- matrix(NA, length(mod1), length(mod1))
for(i in 1:length(mod1)){
  for(j in 1:length(mod1)){
    tempi<- which(mat1$kinase == mod1[i])
    tempj<- which(mat1$kinase2 == mod1[j])
    if(length(intersect(tempi, tempj)) != 0){
      cor1[i,j]<- mat1$cor_cancer[intersect(tempi, tempj)]
    }
  }
} 

fdr1<- matrix(NA, length(mod1), length(mod1))
for(i in 1:length(mod1)){
  for(j in 1:length(mod1)){
    tempi<- which(mat1$kinase == mod1[i])
    tempj<- which(mat1$kinase2 == mod1[j])
    if(length(intersect(tempi, tempj)) != 0){
      fdr1[i,j]<- mat1$eFDR_cancer[intersect(tempi, tempj)]
    }
  }
} 

rownames(cor1)<- colnames(cor1)<- rownames(fdr1)<- colnames(fdr1)<- mod1

ht_figures5e <- Heatmap(as.matrix(cor1), col = colorRamp2(c(-1, 0, 1), c("blue","white", "red")), show_row_names = T,
        show_column_names = T, rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, fill){
          if(!is.na(cor1[i, j]) && (fdr1[i,j] < 0.1)) {
            grid.points(x, y, pch = 16, size = unit(2, "mm"),gp=gpar(fill="black",col="black"))
          }
        }, height = nrow(cor1)*unit(5, "mm"),
        width = ncol(cor1)*unit(5, "mm"),
        column_order = colnames(cor1), #[c(1,3,2,6,7,4,5)]
        row_order = rownames(cor1), #[c(1,3,2,6,7,4,5)]
        cluster_rows = F,
        show_column_dend = F,
        show_row_dend = FALSE, column_title="",
        column_title_gp=grid::gpar(fontsize=12), row_names_side = "left", show_heatmap_legend = T)
save_heatmap_pdf(ht_figures5e, "FigureS5E_module_II_kinase_correlation_heatmap.pdf", width = 5.5, height = 5.5)


############### Figure 5F ###############
score<- read.csv("data/kinase_csn_mat_v2.csv")
rownames(score)<- score$X
score<- score[,-1]
colnames(score)[1:81]<- substr(colnames(score)[1:81],2,nchar(colnames(score)[1:81]))
colnames(score)<- gsub('\\.', '-', colnames(score))

mod1<- score[c(2,3,4,12,21,25,27), colnames(score) %in% colnames(data2.sub)]
mod2<- score[c(5,7,13,16,17,28), colnames(score) %in% colnames(data2.sub)]

p.mod1<- matrix(NA, 7,128)
for(i in 1:7){
  for(j in 1:128){
    p.mod1[i,j]<- pnorm(mod1[i,j], lower.tail = F)
  }
}

p.mod2<- matrix(NA, 6,128)
for(i in 1:6){
  for(j in 1:128){
    p.mod2[i,j]<- pnorm(mod2[i,j], lower.tail = F)
  }
}

colnames(p.mod1)<- colnames(p.mod2)<- colnames(mod1)

mod1.score.p<- mod1.score.z<- mod1.score.combp<- NULL
mod2.score.p<- mod2.score.z<- mod2.score.combp<- NULL
for(k in 1:128){
  tmp1<- sumlog(p.mod1[complete.cases(p.mod1[,k]),k])
  tmp1z<- sumz(p.mod1[complete.cases(p.mod1[,k]),k])
  
  tmp2<- sumlog(p.mod2[complete.cases(p.mod2[,k]),k])
  tmp2z<- sumz(p.mod2[complete.cases(p.mod2[,k]),k])
  
  mod1.score.combp<- c(mod1.score.combp, tmp1$p)
  mod2.score.combp<- c(mod2.score.combp, tmp2$p)
  
  mod1.score.z<- c(mod1.score.z, tmp1z$z)
  mod2.score.z<- c(mod2.score.z, tmp2z$z)
  
  mod1.score.p<- c(mod1.score.p, tmp1z$p)
  mod2.score.p<- c(mod2.score.p, tmp2z$p)
}

df<- rbind(mod1,mod2)

colnames(df)<- colnames(mod1)

df.combp<- rbind(mod1.score.combp,mod2.score.combp)
df.p<- rbind(mod1.score.p,mod2.score.p)

p.mod1.mean<- apply(p.mod1, 2, function(x) mean(x < 0.05))

cell<- c("7316-1746", "7316-1763", "7316-1769", "7316-195","7316-2151","7316-2176","7316-24","7316-3058", "7316-388","7316-445", "7316-85", "7316-913")

cov<- var.tab
colnames(cov)[1]<- "samples"

anno<- cov[cov$samples %in% colnames(mod1),]
anno.ord<- anno[order(match(anno$samples, colnames(mod1))),]
identical(colnames(mod1), anno.ord$samples) # TRUE
anno.ord$cell<- rep(0,128)
anno.ord$cell[which(anno.ord$samples %in% cell)]<- 1
anno.ord$cell<- as.factor(anno.ord$cell)
anno.ord$gender<- as.factor(anno.ord$gender)

p.mod1.mean[is.na(p.mod1.mean)]<- 0
s<- sort(p.mod1.mean,ind=T,decreasing=T)

prop<- p.mod1.mean[s$ix]
mod11<- mod1[,s$ix]
anno.ord1<- anno.ord[s$ix,]

ha <- HeatmapAnnotation(age = anno.ord1$age.seg, sex = anno.ord1$gender, "cell-line"=anno.ord1$cell, "prop"=prop)

ht_figure5f <- Heatmap(as.matrix(mod11), show_row_names = T, col = colorRamp2(c(-1, 0, 5), c("blue","white", "red")),
        show_column_names = F, rect_gp = gpar(col = "white", lwd = 1),
        height = nrow(mod11)*unit(4, "mm"),
        width = ncol(mod11)*unit(2, "mm"),
        column_order = colnames(mod11),
        top_annotation = ha,
        cluster_rows = T,
        cluster_columns = F,
        show_column_dend = F,
        show_row_dend = FALSE, column_title="",
        column_title_gp=grid::gpar(fontsize=12), show_heatmap_legend = T)
save_heatmap_pdf(ht_figure5f, "Figure5F_module_I_kinase_interaction_heatmap.pdf", width = 9, height = 4)

################## Figure 5I #################

drug.data<- read_excel("data/drug_kinase.xlsx")
load("data/data_kea3.RData")

data.kin<- data2.sub[, which(colnames(data2.sub) %in% c("7316-913", "7316-1763", "7316-195", "7316-3058", "7316-388"))]

data.kin<- data.kin[,c(5,1:4)]

drug.data<- drug.data[,c(1,2,6:10)]

plk2.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="PLK2"), which(drug.data$Drug=="BI 2536")),]
atr.drug<- drug.data[which(drug.data$`HOPE target`=="ATR"),]
flt1.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="FLT1"), which(drug.data$Drug=="Foretinib (GSK1363089)")),]
flt1.drug<- flt1.drug[2,]
atm.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="ATM"), which(drug.data$Drug=="VE-822")),]
cdk8.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="CDK8"), which(drug.data$Drug=="Flavopiridol HCl")),]
cdk8.drug<- cdk8.drug[1,]
raf1.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="RAF1"), which(drug.data$Drug=="Vemurafenib (PLX4032- RG7204)")),]
lck.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="LCK"), which(drug.data$Drug=="Ponatinib (AP24534)")),]
lck.drug<- lck.drug[3,]
map1.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="MAP2K1"), which(drug.data$Drug=="Trametinib (GSK1120212)")),]
map1.drug<- map1.drug[3,]
map2.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="MAP2K1"), which(drug.data$Drug=="U0126-EtOH")),]
map3.drug<- drug.data[intersect(which(drug.data$`HOPE target`=="MAP2K1"), which(drug.data$Drug=="MEK162 (ARRY-162_ ARRY-438162)")),]

crisp<- c(t(atm.drug[1,3:7]), t(cdk8.drug[1,3:7]), t(lck.drug[1,3:7]))
kin<- (data.kin[which(rownames(data.kin) %in% c("ATM", "CDK8", "LCK")),])
df0<- cbind.data.frame("crisp"=c((crisp)), "kin"=c(t(kin))) 
df0$kinase<- rep(c("ATM", "CDK8", "LCK"), each=5)

lm_fit <- lm(df0$crisp ~ df0$kin, data=df0)
predicted_df <- data.frame("pred" = predict(lm_fit, df0), "kinase"=df0$kin)
cor.test(df0$crisp,df0$kin, method = "spearman")

p_figure5i_overall <- ggplot(df0, aes(kin,crisp, col=kinase)) +
  geom_point(aes(col=kinase)) +
  scale_colour_manual(values = c("blue", "red", "green")) +
  xlab("Kinase activity scores") +
  ylab("Drug data") +
  geom_line(color='black',data = predicted_df, aes(y=pred, x=kinase), size=1) +
  theme_pubr() +
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 

################ Figure 5I ###############
crisp<- c(t(atm.drug[1,3:7]), t(cdk8.drug[1,3:7]), t(lck.drug[1,3:7]))
kin<- (data.kin[which(rownames(data.kin) %in% c("ATM", "CDK8", "LCK")),])
df0<- cbind.data.frame("crisp"=c((crisp)), "kin"=c(t(kin))) 
df0$kinase<- rep(c("ATM", "CDK8", "LCK"), each=5)

lm_fit0 <- lm(df0$crisp[df0$kinase=="ATM"] ~ df0$kin[df0$kinase=="ATM"], data=df0)
predicted_df0 <- data.frame("pred" = predict(lm_fit0, df0[df0$kinase=="ATM",]), "ATM"=df0$kin[df0$kinase=="ATM"])
atm.cr<- cor.test(df0$crisp[df0$kinase=="ATM"],df0$kin[df0$kinase=="ATM"], method = "spearman")

lm_fit1 <- lm(df0$crisp[df0$kinase=="CDK8"] ~ df0$kin[df0$kinase=="CDK8"], data=df0)
predicted_df1 <- data.frame("pred" = predict(lm_fit1, df0[df0$kinase=="CDK8",]), "CDK8"=df0$kin[df0$kinase=="CDK8"])
cdk8.cr<- cor.test(df0$crisp[df0$kinase=="CDK8"],df0$kin[df0$kinase=="CDK8"], method = "pearson")

lm_fit2 <- lm(df0$crisp[df0$kinase=="LCK"] ~ df0$kin[df0$kinase=="LCK"], data=df0)
predicted_df2 <- data.frame("pred" = predict(lm_fit2, df0[df0$kinase=="LCK",]), "LCK"=df0$kin[df0$kinase=="LCK"])
lck.cr<- cor.test(df0$crisp[df0$kinase=="LCK"],df0$kin[df0$kinase=="LCK"], method = "spearman")

p_figure5i <- ggplot(df0, aes(kin,crisp, col=kinase)) + geom_point(aes(col=kinase)) + scale_colour_manual(values = c("blue", "red", "green")) + xlab("Kinase activity scores") + ylab("Drug data") + geom_line(color='blue',data = predicted_df0, aes(y=pred, x=ATM), size=1)+ geom_line(color='red',data = predicted_df1, aes(y=pred, x=CDK8), size=1) + geom_line(color='green',data = predicted_df2, aes(y=pred, x=LCK), size=1) + theme_pubr() + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 
save_ggplot_pdf(p_figure5i, "Figure5I_ATM_CDK8_LCK_drug_response_scatter.pdf", width = 6, height = 5)

################ Figure S5J ###############
crisp<- c(t(flt1.drug[1,3:7]), t(map2.drug[1,3:7]), t(raf1.drug[1,3:7]))
kin<- (data.kin[which(rownames(data.kin) %in% c("FLT1", "MAP2K1", "RAF1")),])
df0<- cbind.data.frame("crisp"=c((crisp)), "kin"=c(t(kin))) 
df0$kinase<- rep(c("FLT1", "MAP2K1", "RAF1"), each=5)

lm_fit0 <- lm(df0$crisp[df0$kinase=="FLT1"] ~ df0$kin[df0$kinase=="FLT1"], data=df0)
predicted_df0 <- data.frame("pred" = predict(lm_fit0, df0[df0$kinase=="FLT1",]), "FLT1"=df0$kin[df0$kinase=="FLT1"])
flt1.cr<- cor.test(df0$crisp[df0$kinase=="FLT1"],df0$kin[df0$kinase=="FLT1"], method = "spearman")

lm_fit1 <- lm(df0$crisp[df0$kinase=="MAP2K1"] ~ df0$kin[df0$kinase=="MAP2K1"], data=df0)
predicted_df1 <- data.frame("pred" = predict(lm_fit1, df0[df0$kinase=="MAP2K1",]), "MAP2K1"=df0$kin[df0$kinase=="MAP2K1"])
map2.cr<- cor.test(df0$crisp[df0$kinase=="MAP2K1"],df0$kin[df0$kinase=="MAP2K1"], method = "spearman")

lm_fit2 <- lm(df0$crisp[df0$kinase=="RAF1"] ~ df0$kin[df0$kinase=="RAF1"], data=df0)
predicted_df2 <- data.frame("pred" = predict(lm_fit2, df0[df0$kinase=="RAF1",]), "RAF1"=df0$kin[df0$kinase=="RAF1"])
raf1.cr<- cor.test(df0$crisp[df0$kinase=="RAF1"],df0$kin[df0$kinase=="RAF1"], method = "spearman")

p_figures5j <- ggplot(df0, aes(kin,crisp, col=kinase)) + geom_point(aes(col=kinase)) + scale_colour_manual(values = c("blue", "red", "green")) + xlab("Kinase activity scores") + ylab("Drug data") + geom_line(color='blue',data = predicted_df0, aes(y=pred, x=FLT1), size=1)+ geom_line(color='red',data = predicted_df1, aes(y=pred, x=MAP2K1), size=1) + geom_line(color='green',data = predicted_df2, aes(y=pred, x=RAF1), size=1) + theme_pubr() + theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),legend.title = element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15), legend.direction = "vertical", legend.position = "right") 
save_ggplot_pdf(p_figures5j, "FigureS5J_FLT1_MAP2K1_RAF1_drug_response_scatter.pdf", width = 6, height = 5)

################## For Figure 5H, S5I #################
crispr.data<- read_excel("data/crispr_kinase.xlsx")

mat<- as.matrix(crispr.data[,5:8])
colnames(mat)<- gsub("_", "-", colnames(mat))

rownames(mat)<- crispr.data$`HOPE target - Zscore`

data.kin<- data2.sub[rownames(data2.sub) %in% c("ATM", "CDK8", "LCK", "FLT1", "MAP2K1", "RAF1"), which(colnames(data2.sub) %in% c("7316-1746", "7316-195", "7316-3058", "7316-388"))]

mat.sub<- mat[which(rownames(mat) %in% rownames(data.kin)),]
mat.sub.ord<- mat.sub[order(match(rownames(mat.sub), rownames(data.kin))),]  
identical(rownames(data.kin), rownames(mat.sub.ord))
identical(colnames(data.kin), colnames(mat.sub.ord))

cor.vec<- NULL
p.vec<- NULL
for(i in 1:6){
  cr<- cor.test(data.kin[i,], mat.sub.ord[i,], method="spearman")
  cor.vec<- c(cor.vec, cr$estimate)
  p.vec<- c(p.vec, cr$p.value)
}

tab<- cbind.data.frame("Kinases"=rownames(data.kin), "cor.crispr"=cor.vec)
dat1<- cbind.data.frame(tab[c(1,2,4),], "cor.drug"=c(atm.cr$estimate, cdk8.cr$estimate, lck.cr$estimate))
dat2<- cbind.data.frame(tab[-c(1,2,4),], "cor.drug"=c(flt1.cr$estimate, map2.cr$estimate, raf1.cr$estimate))

#%%%%%%%%%% Figure 5H %%%%%%%%%%
df1<- cbind.data.frame("kinase" = rep(c("ATM", "CDK8", "LCK"), 2), "cor"=c(-dat1$cor.crispr, -dat1$cor.drug))
df1$type<- rep(c("CRISPR", "Drug_inhibition"), each=3)
p_figure5h <- ggplot(df1, aes(fill=type, y=cor, x=kinase)) + 
  geom_bar(position="dodge", stat="identity")
save_ggplot_pdf(p_figure5h, "Figure5H_ATM_CDK8_LCK_validation_barplot.pdf", width = 5, height = 4)
# %%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%% Figure S5I %%%%%%%%%%
df2<- cbind.data.frame("kinase" = rep(c("FLT1", "MAP2K1", "RAF1"), 2), "cor"=c(-dat2$cor.crispr, -dat2$cor.drug))
df2$type<- rep(c("CRISPR", "Drug_inhibition"), each=3)
p_figures5i <- ggplot(df2, aes(fill=type, y=cor, x=kinase)) + 
  geom_bar(position="dodge", stat="identity")
save_ggplot_pdf(p_figures5i, "FigureS5I_FLT1_MAP2K1_RAF1_validation_barplot.pdf", width = 5, height = 4)
#%%%%%%%%%%%%%%%%%%%%%%
###### Figure 5E ##########
load("data/cov_kinase.RData")
load("data/data_kea3.RData")
load("data/data_protein_trans_mut.RData")

prot<- data4.sub[,colnames(data4.sub) %in% colnames(data2.sub)]

identical(colnames(prot), colnames(data2.sub))
identical(var.tab$sample_id, colnames(data2.sub)) # TRUE

load("data/pathway_all_db_11292023.rdata")
go<- c(pathway.all.db$HALLMARK, pathway.all.db$KEGG, pathway.all.db$REACTOME, pathway.all.db$BIOCARTA, pathway.all.db$GOBP, pathway.all.db$GOMF, pathway.all.db$MITO3) 

oxphos<- prot[rownames(prot) %in% go$HALLMARK_OXIDATIVE_PHOSPHORYLATION, var.tab$IDH1==0]
var.tab1<- var.tab[var.tab$IDH1==0,]

identical(var.tab1$sample_id, colnames(oxphos)) # TRUE

var.tab1$location<- as.factor(var.tab1$location)
var.tab1$grade<- as.factor(var.tab1$grade)
var.tab1$age.seg<- as.factor(var.tab1$age.seg)

b.m<- matrix(NA,nrow=nrow(oxphos), ncol=1)
p.m<- matrix(NA,nrow=nrow(oxphos), ncol=1)
for(i in 1:dim(oxphos)[1]){
  cur <- as.numeric(oxphos[i,])
  cdk8<- as.numeric(data2.sub[rownames(data2.sub)=="CDK8",var.tab$IDH1==0])
  idh1<- as.numeric(oxphos[rownames(oxphos)=="IDH1",])
  
  temp<- lm(cur ~ cdk8 + age.seg + gender + location + grade, data = var.tab1)

  b.m[i,]<- summary(temp)$coef[2,1]
  p.m[i,]<- summary(temp)$coef[2,4]
}

fdr1<- p.adjust(p.m[,1], method = "BH", n=dim(p.m)[1])

dat<- cbind.data.frame("genes" = rownames(oxphos), "logp" = -log(p.m[,1],10), "FC" = b.m[,1])
#dat<- cbind.data.frame("genes" = tab$genes, "logp" = -log(tab$cdk8.p,10), "FC"=tab$cdk8.coef)
sum(dat$logp > 1)
dat$genes.text<- rep("", dim(dat)[1])
dat$genes.text[dat$logp > 1]<- dat$genes[dat$logp > 1]

p_figures5h <- ggplot(dat, aes(x=FC, y=logp, label=genes.text)) +
  geom_vline(xintercept = 0, linetype = 2, col="black") + geom_hline(yintercept = 1, linetype = 2,col="black") + geom_point(col="grey") + geom_text_repel(size=4,col="blue", max.overlaps = 108) + 
  theme_classic() +
  labs(x="Fold Change (CDK8)",y="-log_10 (pvalue)")  +
  theme(axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),  legend.text=element_text(colour="black",size=18), legend.title=element_text(colour="black",size=20)) 
save_ggplot_pdf(p_figures5h, "Figure5E_CDK8_oxphos_volcano.pdf", width = 6, height = 5)

###### Figure S5G ##########
kids<- read.csv("data/kids_first.tsv", sep="\t")
clinical<- read.csv("data/clinical_table.txt", sep="\t")
data1<- kids[,-c(1:1)]
data1<- as.matrix(data1)
rownames(data1)<- kids$ApprovedGeneSymbol
colnames(data1)<- substr(colnames(data1),2,nchar(colnames(data1)))
colnames(data1)<- gsub('\\.', '-', colnames(data1))
data1<- data1[,-which(colnames(data1) %in% c("7316-1763", "7316-1769", "7316-2146", "7316-942"))]
clinical.sub<- clinical[clinical$sample_id %in% colnames(data1),]
clinical.sub.ord<- clinical.sub[order(match(clinical.sub$sample_id, colnames(data1))),]
identical(clinical.sub.ord$sample_id, colnames(data1))

df<- cbind.data.frame("mut"=clinical.sub.ord$molecular_subtype,"gender"=clinical.sub.ord$reported_gender, "age"=clinical.sub.ord$age_at_diagnosis_days/365, "tumor"=clinical.sub.ord$tumor_descriptor, "OS_days"=(clinical.sub.ord$age_last_update_days - clinical.sub.ord$age_at_diagnosis_days), "OS_days.old"=clinical.sub.ord$OS_days, "status"=clinical.sub.ord$OS_status, t(data1))

df$status1<- ifelse(df$status=="DECEASED",1,0)

ind<- which(df$OS_days > 2000)
df$OS_days[ind]<- 2000
df$status1[ind]<- 0

df$age.group<- rep(NA, dim(df)[1])
df$age.group[intersect(which(df$age >= 0), which(df$age < 15))]<- "PED"
df$age.group[intersect(which(df$age >= 15), which(df$age < 26))]<- "ADOL"

df$mut1<- rep(0, dim(df)[1])
df$mut1[which(grepl("H3 K28", df$mut)==T)]<- 1

df$tum <- rep(0, dim(df)[1])
df$tum[which(grepl("Initial", df$tumor)==T)]<- 1 
df$tum[which(grepl("Progressive", df$tumor)==T)]<- 1 

df$loc <- clinical.sub.ord$CNS_region
df$loc1<- ifelse(df$loc == "Midline", 0,1)

df$OS_days[which(df$OS_days==0)]<- NA

ind.ado<- which(df$age.group=="ADOL")

df0<- df[,c(2,5,10783:10786)]
df0$cl<- ifelse(df$CDK8 > median(df$CDK8[ind.ado]), "high", "low")
df0$cl<- factor(df0$cl, levels = c("low", "high"))
df0$gender<- as.factor(df0$gender)
df0$age.group<- as.factor(df0$age.group)
df0$mut1<- as.factor(df0$mut1)
df0$tum<- as.factor(df0$tum)
df0.m<- df0[df0$gender=="Male",]

cox.fit.new<-tryCatch(
  coxph(
    Surv(OS_days, status1) ~ cl + mut1 + tum + age.group + gender,
    data = df0,
    control = coxph.control(
      eps = 1e-09,
      toler.chol = .Machine$double.eps ^ 0.75,
      iter.max = 10000000,
      toler.inf = sqrt(1e-09),
      outer.max = 100,
      timefix = TRUE
    )
  ),warning=function(w) {NULL},error=function(e) {NULL})

print(summary(cox.fit.new)$coef)

new.data=data.frame(cl=factor(c("high","low"),levels=c("low","high")))
new.data$mut1=0; new.data$mut1<- as.factor(new.data$mut1)
new.data$tum=0; new.data$tum<- as.factor(new.data$tum)
new.data$age.group="ADOL"
new.data$gender="Female"
row.names(new.data)=c("high protein","low protein")

myfit <- survfit(cox.fit.new, newdata = new.data[,c("cl", "mut1", "tum", "age.group","gender")])

pdf(file = file.path(output_dir, "FigureS5G_CDK8_survival_curve.pdf"), width = 6, height=6, useDingbats=T)
ggsurvplot(myfit, data = new.data, size=.5,
           title="KM Curve for CDK8",
           conf.int = FALSE,
           palette=c("high protein"="red","low protein"="black"),
           risk.table = F,
           ggtheme=theme_classic(),
           gg= theme_classic())
dev.off()
