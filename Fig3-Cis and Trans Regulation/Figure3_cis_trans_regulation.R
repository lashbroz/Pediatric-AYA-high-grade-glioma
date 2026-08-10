#!/usr/bin/env Rscript

# ============================================================
# Figure 3 / Figure S3 | Cis and Trans Regulation
# File: Figure3_cis_trans_regulation.R
#
# Description:
#   Generates mutation, CNV, and cis-regulation panels for the
#   Figure 3 and Figure S3 analyses using public manuscript data
#   and supplementary tables.
#
# Inputs:
#   - ../data/STable1.xlsx
#   - ../data/STable3.xlsx
#   - public molecular matrices in ../data/
#
# Outputs:
#   - Figure3A_ATRX_TP53_mutation_protein_RNA_boxplots.pdf
#   - Figure3F_CNV_protein_boxplots.pdf
#   - FigureS3A_NF1_mutation_protein_RNA_boxplots.pdf
#   - FigureS3B_HGG_NOS_mutation_protein_RNA_boxplots.pdf
#   - FigureS3E_CNV_cis_regulation_circle.pdf
#   - FigureS3F_CNV_RNA_boxplots.pdf
#
# Author: Weiping Ma
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

library(openxlsx)

normalize_sample_ids <- function(x) {
  x <- as.character(x)
  gsub("\\.", "-", gsub("^((X|A\\.|G\\.|P\\.))", "", x))
}

read_public_tsv <- function(filename) {
  read.delim(
    file.path(data_dir, filename),
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "", "NaN")
  )
}

split_annotation_matrix <- function(data, annotation_cols, row_id = NULL) {
  annotation <- data[, annotation_cols, drop = FALSE]
  matrix_data <- data[
    ,
    setdiff(seq_along(data), match(names(annotation), names(data))),
    drop = FALSE
  ]
  matrix_data <- as.data.frame(lapply(matrix_data, function(x) suppressWarnings(as.numeric(x))))
  matrix_data <- as.matrix(matrix_data)

  if (!is.null(row_id)) {
    rownames(matrix_data) <- annotation[[row_id]]
  }

  colnames(matrix_data) <- normalize_sample_ids(colnames(matrix_data))
  list(annotation = annotation, matrix = matrix_data)
}

# read data ################################
data_dir = '../data'
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
script_dir <- resolve_script_dir("Figure3_cis_trans_regulation.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

clinical = read.xlsx(file.path(data_dir, "STable1.xlsx"), sheet = "ClinicalTable")
clinical$id = normalize_sample_ids(clinical$id)
if ("sample_id" %in% colnames(clinical)) {
  clinical$sample_id = normalize_sample_ids(clinical$sample_id)
}

mut.raw = read_public_tsv("cDisc_mutation_10192023.tsv")
mut.split = split_annotation_matrix(
  mut.raw,
  annotation_cols = 1,
  row_id = 'ApprovedGeneSymbol'
)
mut = mut.split$matrix

rna.raw = read_public_tsv("cDisc_rna_coding_10192023.tsv")
rna.split = split_annotation_matrix(
  rna.raw,
  annotation_cols = 1:3,
  row_id = 'ApprovedGeneSymbol'
)
rna = rna.split$matrix

cnv.raw = read_public_tsv("cDisc_CNV_coding_10252023.tsv")
cnv.split = split_annotation_matrix(
  cnv.raw,
  annotation_cols = 1,
  row_id = 'ApprovedGeneSymbol'
)
cnv = cnv.split$matrix
g.cnv = cnv.split$annotation$ApprovedGeneSymbol

protein.raw = read_public_tsv("cDisc_proteome_imputed_data_09152023.tsv")
protein.split = split_annotation_matrix(
  protein.raw,
  annotation_cols = 1:4,
  row_id = 'ApprovedGeneSymbol'
)
protein = protein.split$matrix
g.protein = protein.split$annotation$ApprovedGeneSymbol

## read cnv cis results from supp table 3 ##########

cis.res = read.xlsx('../data/STable3.xlsx',sheet = "cis_cnv_methylation")
colnames(cis.res)[4:8] = paste(colnames(cis.res)[4],cis.res[1,4:8])
colnames(cis.res)[4:8+5] = paste(colnames(cis.res)[4+5],cis.res[1,4:8+5])
colnames(cis.res)[4:8+5*2] = paste(colnames(cis.res)[4+5*2],cis.res[1,4:8+5*2])
colnames(cis.res)[4:8+5*3] = paste(colnames(cis.res)[4+5*3],cis.res[1,4:8+5*3])
cis.res = cis.res[which(cis.res$predictor=='cnv'),]

## read gene location from supp table 1 ############
gene_loc = read.xlsx('../data/STable1.xlsx',sheet = "Gene_Location_Annotation")
# figure 3A #######################

library(ggplot2)
library(ggpubr)

save_pdf_plot <- function(plot, filename, width, height) {
  grDevices::pdf(file.path(output_dir, filename), width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
}

save_pdf_grid <- function(grob, filename, width, height) {
  grDevices::pdf(file.path(output_dir, filename), width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()
}

ss = intersect(intersect(colnames(mut),colnames(rna)),colnames(protein))[-1]
ss = intersect(ss, clinical$id[clinical$cDisc_age<62])
cov = clinical[match(ss,clinical$id),]

df.mut.plot = NULL
for( g in c('ATRX','TP53','NF1'))
{
  # g = 'NF1'
  df.g = data.frame(gene = g,
                    sample = ss,
                    rna = unlist(rna[rownames(rna)==g,ss]),
                    pro = unlist(protein[g.protein==g,ss]),
                    protein.zscore = scale(unlist(protein[g.protein==g,ss])),
                    mut = ifelse(unlist(mut[rownames(mut)==g,ss])==1,yes = 'Mut',no = 'WT'),
                    sex = cov$cDisc_Gender,
                    HOPE = cov$cohort,
                    age = cov$cDisc_age,
                    age.cl = (cov$cDisc_age>15)+(cov$cDisc_age>26)+(cov$cDisc_age>40),
                    ped = cov$cDisc_age<15,
                    ado = (cov$cDisc_age>=15)*(cov$cDisc_age<26)
  )
  
  df.mut.plot = rbind(df.mut.plot,df.g)
}

gg.mut.pro=
  ggplot(df.mut.plot[df.mut.plot$gene!='NF1',],aes(x = mut, y = protein.zscore, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + ylim(-5,5)+
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(gene),
             scales="free", space="free")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",label.y = c(3.5),label.x = c(1.5),method = "wilcox.test")

gg.mut.rna=
  ggplot(df.mut.plot[df.mut.plot$gene!='NF1',],aes(x = mut, y = rna, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + ylim(2,8)+
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(gene),
             scales="free", space="free")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",label.y = c(7),label.x = c(1.5),method = "wilcox.test")

save_pdf_plot(
  ggpubr::ggarrange(
    gg.mut.pro + theme(plot.margin = margin(t = 28, r = 6, b = 6, l = 6)),
    gg.mut.rna + theme(plot.margin = margin(t = 28, r = 6, b = 6, l = 6)),
    ncol = 2,
    labels = c("Protein", "RNA"),
    font.label = list(size = 12, face = "bold")
  ),
  "Figure3A_ATRX_TP53_mutation_protein_RNA_boxplots.pdf",
  width = 7,
  height = 5.9
)

# figure S3A #######################

gg.mut.pro.S=
  ggplot(df.mut.plot[df.mut.plot$gene=='NF1',],aes(x = mut, y = protein.zscore, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + ylim(-5,5)+
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(gene),
             scales="free", space="free")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",label.y = c(3.5),label.x = c(1.5),method = "wilcox.test")

gg.mut.rna.S=
  ggplot(df.mut.plot[df.mut.plot$gene=='NF1',],aes(x = mut, y = rna, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + ylim(2,8)+
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(gene),
             scales="free", space="free")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",label.y = c(7),label.x = c(1.5),method = "wilcox.test")

save_pdf_plot(
  ggpubr::ggarrange(
    gg.mut.pro.S + theme(plot.margin = margin(t = 28, r = 6, b = 6, l = 6)),
    gg.mut.rna.S + theme(plot.margin = margin(t = 28, r = 6, b = 6, l = 6)),
    ncol = 2,
    labels = c("Protein", "RNA"),
    font.label = list(size = 12, face = "bold")
  ),
  "FigureS3A_NF1_mutation_protein_RNA_boxplots.pdf",
  width = 7,
  height = 4.0
)

# figure S3B #######################

s.nos = cov$id[which(cov$Disc_Cancer_Group=='(HGG) High Grade Glioma (not otherwise specified)')]

df.mut.plot.sub = dplyr::filter(df.mut.plot,sample%in%s.nos)

df.mut.plot.sub = rbind(data.frame(gene = df.mut.plot.sub$gene,mut = df.mut.plot.sub$mut,
                                   # cnv.cat = df.mut.plot.sub$
                                   value = df.mut.plot.sub$protein.zscore,type = 'protein.z'),
                        data.frame(gene = df.mut.plot.sub$gene,mut = df.mut.plot.sub$mut,
                                   value = df.mut.plot.sub$rna,type = 'RNA'))



gg.mut.ATRX = 
  ggplot(dplyr::filter(df.mut.plot.sub,gene =='ATRX'),
         aes(x = mut, y = value, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(type),
             scales="free", space="free")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",
                     vjust = 1,
                     # label.y = c(3.5),
                     label.x = c(1.5),method = "wilcox.test")

gg.mut.TP53 = 
  ggplot(dplyr::filter(df.mut.plot.sub,gene =='TP53'),
         aes(x = mut, y = value, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(type),
             scales="free", space="fixed")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",
                     vjust = 1,
                     # label.y = c(3.5),
                     label.x = c(1.5),method = "wilcox.test")

gg.mut.NF1 = 
  ggplot(dplyr::filter(df.mut.plot.sub,gene =='NF1'),
         aes(x = mut, y = value, alpha=0.2,fill=mut)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("red","grey"))+scale_color_manual(values=c("red","grey"))+
  facet_grid(rows = vars(type),
             scales="free", space="fixed")+
  geom_point(position = position_jitterdodge(), alpha=0.8,aes(fill = mut),color = 'black',pch=21)+theme_bw()+
  stat_compare_means(aes(group = mut), label = "p.signif",
                     vjust = 1,
                     # label.y = c(3.5),
                     label.x = c(1.5),method = "wilcox.test")

save_pdf_plot(
  ggpubr::ggarrange(gg.mut.ATRX, gg.mut.TP53, gg.mut.NF1, ncol = 3),
  "FigureS3B_HGG_NOS_mutation_protein_RNA_boxplots.pdf",
  width = 9,
  height = 4.5
)

# figure 3F #######################

library(ggplot2)
library(ggpubr)
library(rstatix)
library(circlize)
library(grid)

### data for plot of GUCY1B1, AUTS2, GOLGA2 
g.temp = c('GUCY1B1', 'AUTS2', 'GOLGA2')

df.cnv.plot.main = NULL

for( g in g.temp)
{
  df.g = data.frame(sample = ss,
                    gene = g,
                    rna = unlist(rna[rownames(rna)==g,ss]),
                    pro = unlist(protein[g.protein==g,ss]),
                    protein.zscore = scale(unlist(protein[g.protein==g,ss])),
                    cnv = unlist(cnv[g.cnv==g,ss]),
                    cnv.cat = (unlist(cnv[g.cnv==g,ss])>sd(cnv[g.cnv==g,ss],na.rm = T))-
                      (unlist(cnv[g.cnv==g,ss])<(-sd(cnv[g.cnv==g,ss],na.rm = T))),
                    sex = cov$cDisc_Gender,
                    HOPE = cov$cohort,
                    age = cov$cDisc_age,
                    age.cl = (cov$cDisc_age>15)+(cov$cDisc_age>26)+(cov$cDisc_age>40),
                    ped = cov$cDisc_age<15,
                    ado = (cov$cDisc_age>=15)*(cov$cDisc_age<26)  
  )
  df.cnv.plot.main = rbind(df.cnv.plot.main,df.g)
}

df.cnv.plot.main$age.cat = factor(c('Ped','Ado','YA','Adult')[df.cnv.plot.main$age.cl+1],levels = c('Ped','Ado','YA','Adult'),ordered = T)
df.cnv.plot.main$age.cat2 = factor(c('Ped','AYA','AYA','Adult')[df.cnv.plot.main$age.cl+1],levels = c('Ped','AYA','Adult'),ordered = T)

df.cnv.plot.main$cnv.category = c('loss','neutral','gain')[df.cnv.plot.main$cnv.cat+2]

df.cnv.plot.main$cnv.category[df.cnv.plot.main$gene=='AUTS2'] = 
  ifelse(df.cnv.plot.main$cnv.category[df.cnv.plot.main$gene=='AUTS2']=='gain',yes = 'gain',no = 'non-gain')

df.cnv.plot.main$cnv.category[df.cnv.plot.main$gene=='GOLGA2'] = 
  ifelse(df.cnv.plot.main$cnv.category[df.cnv.plot.main$gene=='GOLGA2']=='gain',yes = 'gain',no = 'non-gain')

df.cnv.plot.main$cnv.category[df.cnv.plot.main$gene=='GUCY1B1'] = 
  ifelse(df.cnv.plot.main$cnv.category[df.cnv.plot.main$gene=='GUCY1B1']=='loss',yes = 'loss',no = 'non-loss')

df.cnv.plot.main$cnv.category = factor(df.cnv.plot.main$cnv.category,levels = c('loss','non-gain','neutral','non-loss','gain'),ordered = T)

gg.cnv.pro2=
  ggplot(dplyr::filter(df.cnv.plot.main,gene!='AUTS2'),
         aes(x = cnv.category, y = protein.zscore)) + ylim(-4,4)+
  geom_boxplot(aes(fill = cnv.category),width = 0.8,position = position_dodge(preserve = "single")) +
  scale_fill_manual(values=c("purple","grey","grey","orange"))+scale_color_manual(values=c("purple","grey","grey","orange"))+
  # coord_flip()+
  # scale_fill_manual(values=c("red","pink","blue","lightblue"))+
  facet_grid(cols = vars(gene),rows = vars(age.cat2),
             scales="free", space="free")+
  # geom_jitter(size=0.4, alpha=0.9,aes(color = as.character(cnv.cat)))
  geom_point(position = position_jitterdodge(jitter.width = 0.8,dodge.width = 1),
             alpha=0.8,aes(fill = (cnv.category)),color = 'black',pch=21)+theme_bw()+
  stat_compare_means( label = "p.signif")

col_fun_age <- colorRamp2(breaks = 0:5*12,
                          colors = c("#146CF6", "#188AF0", "#00B7D8", "#00D4B0", "#00E54B", "#00F800"))
col.temp = col_fun_age(c(0,30,60))
names(col.temp) = c('Ped','AYA','Adult')

g <- ggplot_gtable(ggplot_build(gg.cnv.pro2))

stripr <- which(grepl('strip-r', g$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col.temp[k]
  k <- k+1
}

save_pdf_grid(
  g,
  "Figure3F_CNV_protein_boxplots.pdf",
  width = 7,
  height = 5
)


# figure S3F #######################

gg.cnv.rna2=
  ggplot(dplyr::filter(df.cnv.plot.main,gene!='AUTS2'),
         aes(x = cnv.category, y = rna)) + ylim(0,7)+
  geom_boxplot(aes(fill = cnv.category),width = 0.8,position = position_dodge(preserve = "single")) +
  scale_fill_manual(values=c("purple","grey","grey","orange"))+scale_color_manual(values=c("purple","grey","grey","orange"))+
  # coord_flip()+
  # scale_fill_manual(values=c("red","pink","blue","lightblue"))+
  facet_grid(cols = vars(gene),rows = vars(age.cat2),
             scales="free", space="free")+
  # geom_jitter(size=0.4, alpha=0.9,aes(color = as.character(cnv.cat)))
  geom_point(position = position_jitterdodge(jitter.width = 0.8,dodge.width = 1),
             alpha=0.8,aes(fill = (cnv.category)),color = 'black',pch=21)+theme_bw()

col_fun_age <- colorRamp2(breaks = 0:5*12,
                          colors = c("#146CF6", "#188AF0", "#00B7D8", "#00D4B0", "#00E54B", "#00F800"))
col.temp = col_fun_age(c(0,30,60))
names(col.temp) = c('Ped','AYA','Adult')

g <- ggplot_gtable(ggplot_build(gg.cnv.rna2))

stripr <- which(grepl('strip-r', g$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col.temp[k]
  k <- k+1
}

save_pdf_grid(
  g,
  "FigureS3F_CNV_RNA_boxplots.pdf",
  width = 7,
  height = 5
)

## t test pvalues of Figure 3F and S3F ################

safe_ttest_p <- function(formula, data) {
  tryCatch(
    t.test(formula, data = data)[[3]],
    error = function(e) NA_real_
  )
}

diagnostic_ttest_pvalues <- data.frame(
  gene = rep(c("GUCY1B1", "GOLGA2", "GOLGA2", "GUCY1B1"), each = 3),
  assay = rep(c("protein", "protein", "RNA", "RNA"), each = 3),
  age.cat2 = rep(c("Ped", "AYA", "Adult"), times = 4),
  p.value = c(
    safe_ttest_p(pro ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GUCY1B1" & age.cat2 == "Ped")),
    safe_ttest_p(pro ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GUCY1B1" & age.cat2 == "AYA")),
    safe_ttest_p(pro ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GUCY1B1" & age.cat2 == "Adult")),
    safe_ttest_p(pro ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GOLGA2" & age.cat2 == "Ped")),
    safe_ttest_p(pro ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GOLGA2" & age.cat2 == "AYA")),
    safe_ttest_p(pro ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GOLGA2" & age.cat2 == "Adult")),
    safe_ttest_p(rna ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GOLGA2" & age.cat2 == "Ped")),
    safe_ttest_p(rna ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GOLGA2" & age.cat2 == "AYA")),
    safe_ttest_p(rna ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GOLGA2" & age.cat2 == "Adult")),
    safe_ttest_p(rna ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GUCY1B1" & age.cat2 == "Ped")),
    safe_ttest_p(rna ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GUCY1B1" & age.cat2 == "AYA")),
    safe_ttest_p(rna ~ cnv.category, data = dplyr::filter(df.cnv.plot.main, gene == "GUCY1B1" & age.cat2 == "Adult"))
  )
)
print(diagnostic_ttest_pvalues)


# figure S3E circle ###################

library(circlize)

cnv = cnv[,colnames(cnv)%in%clinical$id[!clinical$cDisc_age>62]]
ind.cnv.keep = which(rowSums(is.na(cnv[,]))<10)
g.cnv = g.cnv[ind.cnv.keep]
cnv = cnv[ind.cnv.keep,]

# cut.cnv = 0.2
cut.cnv = apply(cnv,1,sd,na.rm = T)

prop_gain = aggregate(t(cnv>cut.cnv), 
                      by = list(age.cl = clinical$cDisc_age_class_name_derived[match(colnames(cnv),clinical$id)]),
                      FUN = mean,na.rm = T)
prop_loss = aggregate(t((-cnv)>cut.cnv), 
                      by = list(age.cl = clinical$cDisc_age_class_name_derived[match(colnames(cnv),clinical$id)]),
                      FUN = mean,na.rm = T)

rownames(prop_gain) = prop_gain[,1]
rownames(prop_loss) = prop_loss[,1]
# prop_gain[is.na(prop_gain)] = 0
# prop_loss[is.na(prop_loss)] = 0

prop_gain = prop_gain[,-1]
prop_loss = prop_loss[,-1]

colnames(prop_gain) = g.cnv
colnames(prop_loss) = g.cnv

gene.loc = gene_loc[,c(1:3,8)]
colnames(gene.loc) = c('chr',     'start',       'end','gene')
gene.loc = gene.loc[match(g.cnv,gene.loc$gene),]

gene.loc$start = as.integer(gene.loc$start)
gene.loc$end = as.integer(gene.loc$end)

# prop_gain.sub = prop_gain[,match(gene.loc$gene,colnames(prop_gain))]
# prop_loss.sub = prop_loss[,match(gene.loc$gene,colnames(prop_loss))]

# mean(gene.loc$gene==colnames(prop_gain))
# mean(gene.loc$gene==colnames(prop_loss))

f.col.age = colorRamp2(breaks = c(1,2.5,4),colors = c('blue','grey80','darkgreen'))

col_fun_age <- colorRamp2(breaks = 0:5*12,
                          colors = c("#146CF6", "#188AF0", "#00B7D8", "#00D4B0", "#00E54B", "#00F800"))

# f.col.age(1:4)
col.subtype = c(col_fun_age(c(0,20,40,60)),'black')
names(col.subtype) = c('PED','ADO','YA','ADULT','All')

factors = paste0("chr", c(1:22,'X'))[23]

# ### cascade overall and age groups specific
# c.r = res20.sig$yName.1[((res20.sig$xType=='cnv')*(res20.sig$yType=='rna'))==1]
# c.pr = res20.sig$yName.2[((res20.sig$xType=='cnv')*(res20.sig$yType=='protein'))==1]
# c.ph = res20.sig$yName.3[((res20.sig$xType=='cnv')*(res20.sig$yType=='phospho'))==1]
# 
# c_p.r = res20.sig$yName.1[((res20.sig$xType=='cnv_p')*(res20.sig$yType=='rna'))==1]
# c_p.pr = res20.sig$yName.2[((res20.sig$xType=='cnv_p')*(res20.sig$yType=='protein'))==1]
# c_p.ph = res20.sig$yName.3[((res20.sig$xType=='cnv_p')*(res20.sig$yType=='phospho'))==1]
# 
# c_a.r = res20.sig$yName.1[((res20.sig$xType=='cnv_a')*(res20.sig$yType=='rna'))==1]
# c_a.pr = res20.sig$yName.2[((res20.sig$xType=='cnv_a')*(res20.sig$yType=='protein'))==1]
# c_a.ph = res20.sig$yName.3[((res20.sig$xType=='cnv_a')*(res20.sig$yType=='phospho'))==1]
# 
# c_y.r = res20.sig$yName.1[((res20.sig$xType=='cnv_y')*(res20.sig$yType=='rna'))==1]
# c_y.pr = res20.sig$yName.2[((res20.sig$xType=='cnv_y')*(res20.sig$yType=='protein'))==1]
# c_y.ph = res20.sig$yName.3[((res20.sig$xType=='cnv_y')*(res20.sig$yType=='phospho'))==1]
# 
# #### protein/RNA cascade 
# ca.c = intersect(c.r,c.pr)
# ca.c_p = intersect(c_p.r,c_p.pr)
# ca.c_a = intersect(c_a.r,c_a.pr)
# ca.c_y= intersect(c_y.r,c_y.pr)
# 
# #### protein/RNA/phospho cascade 
# ca.c3 = intersect(intersect(c.r,c.pr),c.ph)
# ca.c3_p = intersect(intersect(c_p.r,c_p.pr),c_p.ph)
# ca.c3_a = intersect(intersect(c_a.r,c_a.pr),c_a.ph)
# ca.c3_y= intersect(intersect(c_y.r,c_y.pr),c_y.ph)
# 
# setdiff(ca.c,cis.res$gene[cis.res$`All.sample.cis RNA_Protein_Cascade`=='yes'])
# setdiff(cis.res$gene[cis.res$`All.sample.cis RNA_Protein_Cascade`=='yes'],ca.c)
# 
# setdiff(ca.c_p,cis.res$gene[cis.res$`Pediatric.specific.cis RNA_Protein_Cascade_Ped`=='yes'])
# setdiff(cis.res$gene[cis.res$`Pediatric.specific.cis RNA_Protein_Cascade_Ped`=='yes'],ca.c_p)
# 
# setdiff(ca.c_y,cis.res$gene[cis.res$`Young.Adult.specific.cis RNA_Protein_Cascade_YA`=='yes'])
# setdiff(cis.res$gene[cis.res$`Young.Adult.specific.cis RNA_Protein_Cascade_YA`=='yes'],ca.c_y)
# 
# setdiff(ca.c_a,cis.res$gene[cis.res$`Adolescent.specific.cis RNA_Protein_Cascade_Ado`=='yes'])
# setdiff(cis.res$gene[cis.res$`Adolescent.specific.cis RNA_Protein_Cascade_Ado`=='yes'],ca.c_a)
# 
# 
# setdiff(ca.c3,cis.res$gene[cis.res$`All.sample.cis RNA_Protein_Phospho_Cascade`=='yes'])
# setdiff(cis.res$gene[cis.res$`All.sample.cis RNA_Protein_Phospho_Cascade`=='yes'],ca.c3)
# 
# setdiff(ca.c3_p,cis.res$gene[cis.res$`Pediatric.specific.cis RNA_Protein_Phospho_Cascade_Ped`=='yes'])
# setdiff(cis.res$gene[cis.res$`Pediatric.specific.cis RNA_Protein_Phospho_Cascade_Ped`=='yes'],ca.c3_p)
# 
# setdiff(ca.c3_y,cis.res$gene[cis.res$`Young.Adult.specific.cis RNA_Protein_Phospho_Cascade_YA`=='yes'])
# setdiff(cis.res$gene[cis.res$`Young.Adult.specific.cis RNA_Protein_Phospho_Cascade_YA`=='yes'],ca.c3_y)
# 
# setdiff(ca.c3_a,cis.res$gene[cis.res$`Adolescent.specific.cis RNA_Protein_Phospho_Cascade_Ado`=='yes'])
# setdiff(cis.res$gene[cis.res$`Adolescent.specific.cis RNA_Protein_Phospho_Cascade_Ado`=='yes'],ca.c3_a)

ca.c = cis.res$gene[cis.res$`All.sample.cis RNA_Protein_Cascade`=='yes']
ca.c_p = cis.res$gene[cis.res$`Pediatric.specific.cis RNA_Protein_Cascade_Ped`=='yes']
ca.c_y = cis.res$gene[cis.res$`Young.Adult.specific.cis RNA_Protein_Cascade_YA`=='yes']
ca.c_a = cis.res$gene[cis.res$`Adolescent.specific.cis RNA_Protein_Cascade_Ado`=='yes']

ca.c3 = cis.res$gene[cis.res$`All.sample.cis RNA_Protein_Phospho_Cascade`=='yes']
ca.c3_p = cis.res$gene[cis.res$`Pediatric.specific.cis RNA_Protein_Phospho_Cascade_Ped`=='yes']
ca.c3_y = cis.res$gene[cis.res$`Young.Adult.specific.cis RNA_Protein_Phospho_Cascade_YA`=='yes']
ca.c3_a = cis.res$gene[cis.res$`Adolescent.specific.cis RNA_Protein_Phospho_Cascade_Ado`=='yes']

### cnv +; cnv_p +
# g.temp = intersect(ca.c_p,ca.c)[sig.p$est[match(intersect(ca.c_p,ca.c),sig.p$xName)]>1]

### cnv ns; cnv_p -

# sample_p = intersect(colnames(cnv),clinical$id[clinical$cDisc_age_class_name_derived=='PED'])
# sample_a = intersect(colnames(cnv),clinical$id[clinical$cDisc_age_class_name_derived=='ADO'])
# sample_y = intersect(colnames(cnv),clinical$id[clinical$cDisc_age_class_name_derived=='YA'])

up.cna = rowMeans(cnv>cut.cnv,na.rm = T)
down.cna = rowMeans((-cnv)>cut.cnv,na.rm = T)
# up.cna_p = rowMeans(cnv[,sample_p]>cut.cnv,na.rm = T)
# down.cna_p = rowMeans((-cnv[,sample_p])>cut.cnv,na.rm = T)
# up.cna_a = rowMeans(cnv[,sample_a]>cut.cnv,na.rm = T)
# down.cna_a = rowMeans((-cnv[,sample_a])>cut.cnv,na.rm = T)
# up.cna_y = rowMeans(cnv[,sample_y]>cut.cnv,na.rm = T)
# down.cna_y = rowMeans((-cnv[,sample_y])>cut.cnv,na.rm = T)

cut.sample = 0.20
active.c = g.cnv[(up.cna+down.cna)>cut.sample]
active.c_p = g.cnv[(prop_gain+prop_loss)['PED',]>cut.sample]
active.c_a = g.cnv[(prop_gain+prop_loss)['ADO',]>cut.sample]
active.c_y = g.cnv[(prop_gain+prop_loss)['YA',]>cut.sample]

df.KS.pro = rbind(data.frame(gene = intersect(ca.c,active.c), subtype = 4),
                  data.frame(gene = intersect(ca.c_p,active.c_p), subtype = 1),
                  data.frame(gene = intersect(ca.c_a,active.c_a), subtype = 2),
                  data.frame(gene = intersect(ca.c_y,active.c_y), subtype = 3))

df.KS.pho = rbind(data.frame(gene = intersect(ca.c3,active.c), subtype = 4),
                  data.frame(gene = intersect(ca.c3_p,active.c_p), subtype = 1),
                  data.frame(gene = intersect(ca.c3_a,active.c_a), subtype = 2),
                  data.frame(gene = intersect(ca.c3_y,active.c_y), subtype = 3))

df.KS.pro = cbind(gene.loc[match(df.KS.pro$gene,gene.loc$gene),],df.KS.pro)
df.KS.pro = df.KS.pro[!df.KS.pro$chr=='chrX',];
bed.W = list(df.KS.pro[df.KS.pro$subtype==3,c(1:3,6)],df.KS.pro[df.KS.pro$subtype==2,c(1:3,6)],
             df.KS.pro[df.KS.pro$subtype==1,c(1:3,6)],df.KS.pro[df.KS.pro$subtype==4,c(1:3,6)])

df.KS.pho = cbind(gene.loc[match(df.KS.pho$gene,gene.loc$gene),],df.KS.pho)
df.KS.pho = df.KS.pho[!df.KS.pho$chr=='chrX',];
bed.P = list(df.KS.pho[df.KS.pho$subtype==3,c(1:3,6)],df.KS.pho[df.KS.pho$subtype==2,c(1:3,6)],
             df.KS.pho[df.KS.pho$subtype==1,c(1:3,6)],df.KS.pho[df.KS.pho$subtype==4,c(1:3,6)])


k.list = c("ADCK1", "AKT2", "BRAF", "BRSK1", "CAMK1D", "CAMK2G", "CDC42BPB","CDK13", "CDK14", 
"CDK5", "CHUK", "CLK2", "COQ8B", "CSNK1G2","DCLK1", "FYN", "GRK5", "GSK3A", 
"INSR", "JAK2", "LATS1", "LMTK2", "MAP2K2", "MAP2K7", "MAP3K10","MAP3K4", "MAP3K7",
"MAP4K5", "MAPK8", "NEK9", "PAK4", "PAN3", "PIP5K1A","PIP5K1C","PKN1", "PRKD1", 
"PRKG1", "RIPK1", "SLK", "SRPK2", "STK11", "STK32C", "TYK2", "VRK1", "CSNK2A2",
"DAPK1", "MAP2K5", "PRAG1", "PRKAA1", "PRKCI", "RPS6KA4","SIK3", "CSNK2A1","HCK", 
"IRAK4", "SCYL2", "VRK3", "WNK2", "AKT1", "MARK3")

df.KS.pho.sub = df.KS.pho[df.KS.pho$gene%in%k.list,-5]
site.temp = do.call(rbind,strsplit(cis.res$Phosphosite_name[match(df.KS.pho.sub$gene,cis.res$gene)],split = '_'))[,7]
# site.temp = do.call(rbind,strsplit(res20.sig$yName.4[match(df.KS.pho.sub$gene,res20.sig$xName)],split = '_'))[,7]
site.temp[site.temp=='NP'] = 'NLS'

df.KS.pho.sub$gene = paste(df.KS.pho.sub$gene,site.temp,sep = '_')

df.KS.all = rbind(df.KS.pro[df.KS.pro$gene%in%k.list,-5],
                  df.KS.pho.sub)


# col.subtype

bed.PED = rbind(data.frame(gene.loc[,-4],prop = t(prop_gain['PED',]))[t(prop_gain['PED',])>0,],
                data.frame(gene.loc[,-4],prop = t(-prop_loss['PED',]))[t(prop_loss['PED',])>0,])

bed.ADO = rbind(data.frame(gene.loc[,-4],prop = t(prop_gain['ADO',]))[t(prop_gain['ADO',])>0,],
                data.frame(gene.loc[,-4],prop = t(-prop_loss['ADO',]))[t(prop_loss['ADO',])>0,])

bed.YA = rbind(data.frame(gene.loc[,-4],prop = t(prop_gain['YA',]))[t(prop_gain['YA',])>0,],
               data.frame(gene.loc[,-4],prop = t(-prop_loss['YA',]))[t(prop_loss['YA',])>0,])

bed.ADULT = rbind(data.frame(gene.loc[,-4],prop = t(prop_gain['ADULT',]))[t(prop_gain['ADULT',])>0,],
                  data.frame(gene.loc[,-4],prop = t(-prop_loss['ADULT',]))[t(prop_loss['ADULT',])>0,])

# bed.PED = bed.PED[!grepl('X',bed.PED$chr),]

g.KS = do.call(rbind,strsplit(df.KS.all$gene,split = '_'))[,1]

id.s = c(1,4,6,15,26,29,30,37,42,43,49,50,53,60,66,68,84,88,11,14,16,27,40,41,47,55,64,67,73,77,93,95,32,35,72)

# id.s = unique(c(which(rowSums(s.f.pro[g.KS,10:12]!=0)>0),
#                 which(rowSums(s.m.pro[g.KS,10:12]!=0)>0),
#                 which(rowSums(s.f.rna[g.KS,10:12]!=0)>0),
#                 which(rowSums(s.m.rna[g.KS,10:12]!=0)>0)))

df.KS.all$gene[id.s] = paste0('*',df.KS.all$gene[id.s])

col.bar = c('orange','purple')
# col.bar = c('orange','lightblue')
# circos.yaxis(side = 'right', at = c(-0.5,0,0.5),labels = c('50%',0,'50%') , sector.index = 'chrX', track.index = 6 )

lwd.temp = 1

f = colorRamp2(breaks = c(1:4), colors = col.subtype[c(1:3,5)])

pdf(file.path(output_dir, 'FigureS3E_CNV_cis_regulation_circle.pdf'),height = 10,width = 10)
par(mar = c(0,0,0,0))

circos.par(gap.degree = c(rep(1/2,21),3.5,3.5),track.margin = c(0.003,0.003),start.degree=270-8.3-3.5,
           #   "track.height" = 0.05, 
           cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(  plotType = NULL,
                                #   plotType = c("labels"),
                                #   plotType = c("labels", "axis"),
                                chromosome.index = paste0("chr", c(1:22,'X')))

circos.genomicLabels(df.KS.all, labels.column = 4, side = "outside",col = col.subtype[-4][df.KS.all$subtype],
                     connection_height = 0.035,cex = 0.8,labels_height = 0.15)

circos.genomicTrackPlotRegion(bed.P, bg.col = 'grey90',
                              stack = TRUE,track.height = 0.08,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = f(value[[1]]), ...)
                                i = getI(...)
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              })

circos.genomicTrackPlotRegion(bed.W, bg.col = 'grey90',
                              stack = TRUE,track.height = 0.08,
                              track.margin = c(0.01,0.003),
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = f(value[[1]]), ...)
                                i = getI(...)
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              })


circos.genomicTrack(bed.PED, 
                    ylim = c(-0.8,0.8),track.height = 0.05,
                    bg.border = col.subtype['PED'],
                    bg.lwd = lwd.temp,
                    #                     bg.col = 'grey',
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         lwd.border = 1/10,
                                         #                                          border = F,
                                         border = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]),
                                         col = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "black",lwd = 1/4)
                    })

circos.genomicTrack(bed.ADO, 
                    ylim = c(-0.8,0.8),track.height = 0.05,
                    bg.border = col.subtype['ADO'],
                    bg.lwd = lwd.temp,
                    #                     bg.col = 'grey',
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         lwd.border = 1/10,
                                         #                                          border = F,
                                         border = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]),
                                         col = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "black",lwd = 1/4)
                    })

circos.genomicTrack(bed.YA, 
                    ylim = c(-0.8,0.8),track.height = 0.05,
                    bg.border = col.subtype['YA'],
                    bg.lwd = lwd.temp,
                    #                     bg.col = 'grey',
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         lwd.border = 1/10,
                                         #                                          border = F,
                                         border = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]),
                                         col = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "black",lwd = 1/4)
                    })

circos.genomicTrack(bed.ADULT, 
                    ylim = c(-0.8,0.8),track.height = 0.05,
                    bg.border = col.subtype['ADULT'],
                    bg.lwd = lwd.temp,
                    #                     bg.col = 'grey',
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         lwd.border = 1/10,
                                         #                                          border = F,
                                         border = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]),
                                         col = ifelse(value[[1]] > 0, col.bar[1],col.bar[2]), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "black",lwd = 1/4)
                    })


for(tr in 5:8)
{
  circos.yaxis(side = 'left', at = c(-0.5,0,0.5),labels = c('50',0,'50') , sector.index = 'chr1', track.index = tr ,
               labels.cex = 0.04*(16-tr),tick.length = 400000*(14-tr))
}

circos.trackText(x = rep(155*10^6/2,1) , y = rep(3.3,1),
                 labels = 'RNA/Protein/Phospho',track = 3,
                 cex = 0.6, factors = factors, col = "black", font = 1.5, facing = "bending.inside",niceFacing = TRUE)

circos.trackText(x = rep(155*10^6/2,1) , y = rep(1.9,1),
                 labels = 'Cascade CNV',track = 3,
                 cex = 0.6, factors = factors, col = "black", font = 1, facing = "bending.inside",niceFacing = TRUE)

circos.trackText(x = rep(155*10^6/2,1) , y = rep(3.3,1),
                 labels = 'RNA/Protein',track = 4,
                 cex = 0.6, factors = factors, col = "black", font = 1, facing = "bending.inside",niceFacing = TRUE)

circos.trackText(x = rep(155*10^6/2,1) , y = rep(1.9,1),
                 labels = 'Cascade CNV',track = 4,
                 cex = 0.6, factors = factors, col = "black", font = 1, facing = "bending.inside",niceFacing = TRUE)


circos.trackText(x = rep(155*10^6/2,1) , y = rep(0,1),
                 labels = 'PED',track = 5,
                 cex = 0.8, factors = factors, col = 'black', font = 2, facing = "bending.inside",niceFacing = TRUE)

circos.trackText(x = rep(155*10^6/2,1) , y = rep(0,1),
                 labels = 'ADO',track = 6,
                 cex = 0.8, factors = factors, col = 'black', font = 2, facing = "bending.inside",niceFacing = TRUE)

circos.trackText(x = rep(155*10^6/2,1) , y = rep(0,1),
                 labels = 'YA',track = 7,
                 cex = 0.8, factors = factors, col = 'black', font = 2, facing = "bending.inside",niceFacing = TRUE)

circos.trackText(x = rep(155*10^6/2,1) , y = rep(0,1),
                 labels = 'ADULT',track = 8,
                 cex = 0.8, factors = factors, col = 'black', font = 2, facing = "bending.inside",niceFacing = TRUE)


circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  #   circos.rect(xlim[1], 0, xlim[2], 1, col = 'white')
  circos.text(mean(xlim), weighted.mean(ylim,w = c(3,1)), gsub('','',chr)[!grepl('X',chr)], cex = 0.8, col = "black",
              facing = "clockwise", niceFacing = T,font = 2)
}, track.height = 0.06, bg.border = 'white')

# highlight.sector(sector.index = 'chr14', col = "#FFFF0040", track.index = 5:8)

# legend("bottomleft", pch = 1, legend = "This is the legend")
legend('center',legend = c(names(col.subtype),'% of Samples with CN Gain','% of Samples with CN Loss'), 
       cex = 0.8,pch = rep(15,7),
       #        pch = c(rep(124,4),15,15),
       bty = 'n',
       col = c(col.subtype,col.bar), title = "CNV-RNA/Protein/Phospho\nRegulations"
       #        text.col = col.subtype
)
circos.clear()

dev.off();
