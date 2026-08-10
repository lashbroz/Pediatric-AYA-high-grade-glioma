# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
#
# Purpose:
#   Generate the Figure S1D Discovery-cohort copy-number landscape overview
#   using the public CNV matrix and clinical annotations from the
#   repository-level data directory.

library(openxlsx)
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
script_dir <- resolve_script_dir("FigureS1D_CNV_landscape.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# read data ################################
clinical = read.xlsx('../data/STable1.xlsx', sheet = "ClinicalTable")

cnv.raw = read.delim(
  '../data/cDisc_CNV_coding_10252023.tsv',
  sep = '\t',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
g.cnv = cnv.raw$ApprovedGeneSymbol
cnv = as.matrix(cnv.raw[, setdiff(colnames(cnv.raw), 'ApprovedGeneSymbol'), drop = FALSE])
mode(cnv) = 'numeric'
rownames(cnv) = g.cnv

## read gene location from supp table 1 ############
gene_loc = read.xlsx('../data/STable1.xlsx',sheet = "Gene_Location_Annotation")
gene.loc = gene_loc[,c(1:3,8)]
colnames(gene.loc) = c('chr',     'start',       'end','gene')
gene.loc = gene.loc[match(g.cnv,gene.loc$gene),]
gene_loc_cnv = gene.loc
colnames(gene_loc_cnv)[1] = 'Chromosome'

######################################################
### landscape U40

sample.U40 = intersect(colnames(cnv),clinical$id[!clinical$cDisc_age>40])

# cov.U40 = t(cov)[sample.U40,-6]
# cov.U40 = data.frame(cov.U40,Ped = cov.U40[,2]<15,Adult = cov.U40[,2]>26)
# # cov.U40 = data.frame(cov.U40,Ped = cov.U40[,2]<9.5,Adult = cov.U40[,2]>25)
# cov.U40 = cov.U40[,-2]
# cov.U40$cl1 = cl.pro[match(rownames(cov.U40),cl.pro$id),3]==1
# cov.U40$cl2 = cl.pro[match(rownames(cov.U40),cl.pro$id),3]==2
# cov.U40$cl3 = cl.pro[match(rownames(cov.U40),cl.pro$id),3]==3

# cnv.U40 = as.matrix(cnv[,sample.U40])
cnv.U40 = as.matrix(cnv[,sample.U40])
id.valid = which(rowSums(is.na(cnv.U40))<20)

cnv.U40 = cnv.U40[id.valid,]
gene_loc_U40 = gene_loc_cnv[id.valid,]

# arm_loc = read.csv('D:/work/chromosome_arms.csv')
# arm_loc$acen = as.numeric(gsub(',','',arm_loc$acen))
# arm_loc$chr = paste0('chr',arm_loc$chr)
# 
# gene_loc_U40$arm = '_q'
# gene_loc_U40$arm[gene_loc_U40$end<(arm_loc$acen[match(gene_loc_U40$Chromosome,arm_loc$chr)])]   = '_p'
# gene_loc_U40$arm = paste0(gene_loc_U40$Chromosome,gene_loc_U40$arm)

k = 10^ceiling(log10(max(gene_loc_cnv$end,na.rm = T)))
gene_loc_U40$loc = gene_loc_U40$end+k*(as.numeric(gsub('Y','24',gsub('X','23',(gsub('chr','',gene_loc_U40$Chromosome))))))

# head(gene_loc_U40)
loc.na = is.na(gene_loc_U40$loc)
loc.y = grepl('Y',gene_loc_U40$Chromosome)
loc.x = grepl('X',gene_loc_U40$Chromosome)
cnv.U40 = cnv.U40[(loc.na+loc.y+loc.x)==0,]
gene_loc_U40 = gene_loc_U40[(loc.na+loc.y+loc.x)==0,]

# sample.new%in%colnames(cnv.U40)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

col.rdt = brewer.pal(4,'Set3')
names(col.rdt) = 1:4

col_fun_age <- colorRamp2(breaks = 0:5*12,
                          colors = c("#146CF6", "#188AF0", "#00B7D8", "#00D4B0", "#00E54B", "#00F800"))


set.seed(1234)
ha1 = rowAnnotation(df = data.frame(
  # rdt = as.character(4-cov.U40$cl1*3-cov.U40$cl2*2-cov.U40$cl3*1),
  # new.sample = colnames(cnv.U40)%in%sample.new,
  age = clinical$cDisc_age[match(colnames(cnv.U40),clinical$id)]),
  col = list(age = col_fun_age))

# rownames(cnv.U40) = gene_loc_U40$Approved.symbol

set.seed(1234)
h.U40.r = Heatmap(t(cnv.U40),name = 'CNV',
                  col = colorRamp2(breaks = c(-1,0,1),colors = c('purple','white','orange')),
                  cluster_columns = F,
                  show_row_names = F,
                  show_column_names = F,
                  column_order = order(gene_loc_U40$loc),
                  column_split = as.numeric(gsub('chr','',gene_loc_U40$Chromosome)),
                  # column_order = order(gene_loc_U40[which(gene_loc_U40$Chromosome=='chr10'),2]),
                  # row_order = order(apply(cnv.U40[which(gene_loc_U40$Chromosome=='chr10'),],2,mean,na.rm = T)),
                  right_annotation = ha1
                  ,use_raster = T
                  # ,top_annotation = ha2
)

pdf(file.path(output_dir, 'FigureS1D_CNV_landscape.pdf'),height = 6,width = 12)
draw(h.U40.r)
dev.off()
