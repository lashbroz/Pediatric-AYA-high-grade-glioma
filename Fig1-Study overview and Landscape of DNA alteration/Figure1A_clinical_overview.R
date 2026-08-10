#!/usr/bin/env Rscript

# ============================================================
# Figure 1A | Discovery Cohort Clinical Overview
# File: Figure1A_clinical_overview.R
#
# Description:
#   Generates the circular clinical annotation overview for the
#   Discovery cohort using the public clinical table in STable1.
#
# Inputs:
#   - ../data/STable1.xlsx, sheet ClinicalTable
#
# Outputs:
#   - Figure1A_clinical_overview.pdf
#
# Author: Weiping Ma
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

# Load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse) 
  library(ComplexHeatmap)
  library(circlize) 
  library(openxlsx)
})

# Define directories
data_dir <- '../data'
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
script_dir <- resolve_script_dir("Figure1A_clinical_overview.R")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load clinical data
annot <- read.xlsx(file.path(data_dir, "STable1.xlsx"), sheet = "ClinicalTable")

# Filter out missing diagnoses
annot <- annot %>% filter(!is.na(Disc_initial_diagnosis))

# Fix diagnosis type field
annot <- annot %>%
  mutate(cDisc_diagnosis_type = case_when(
    cDisc_diagnosis_type == "Primary" ~ "Initial CNS Tumor",
    TRUE ~ as.character(cDisc_diagnosis_type)
  ))

annot$Diagnosis <- annot$Disc_Diagnosis

# 2. Prepare dataframe for plotting
plot_df <- annot %>%
  dplyr::rename("Age" = "cDisc_age_class_derived",
                "Age_at_Initial_Diagnosis" = "Disc_Age_at_Initial_Diagnosis",
                "Gender" = "cDisc_Gender",
                "Diagnosis Type" = "cDisc_diagnosis_type",
                "Treatment" = "cDisc_treat_status",
                "Tumor Location" = "cDisc_tumor_loc") %>%
  dplyr::mutate(Age_at_Initial_Diagnosis = Age_at_Initial_Diagnosis/365) %>%
  dplyr::select(id, Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, 
                `Diagnosis Type`, Treatment, `Tumor Location`) %>%
  unique() %>%
  column_to_rownames('id') %>%
  dplyr::arrange(Age, Age_at_Initial_Diagnosis, Gender, Diagnosis, 
                 `Diagnosis Type`, Treatment, `Tumor Location`)

# Merge some age groups
plot_df <- plot_df %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Age = ifelse(Age %in% c("(15,26]", "(26,40]"), "(15,40]", Age))

# Splitting variable for circular heatmap
split <- factor(plot_df$Age, levels = c("[0,15]", "(15,40]"))

# 3. Add dummy rows and finalize dataframe
# Add 6 empty slots to the matrix
k <- 6
plot_df2 <- rbind(matrix(NA, ncol = ncol(plot_df), nrow = k, 
                         dimnames = list(1:k, colnames(plot_df))),
                  plot_df)
plot_df2$Age[is.na(plot_df2$Age)] <- '[0]'

split2 <- factor(plot_df2$Age, levels = c('[0]', "[0,15]", "(15,40]"))

# Clean tumor location labels
plot_df2$`Tumor Location`[plot_df2$`Tumor Location` == "Other/Multiple locations/NOS"] <- 'Other/Multiple/NOS'

# Function to add sample sizes
add.sample.size <- function(x){
  tab.n <- table(x)
  for(i in names(tab.n)) {
    x[x == i] <- paste0(i, ' (n=', tab.n[i], ')')
  }
  return(x)
}

# Apply function
plot_df2$Age <- add.sample.size(plot_df2$Age)
plot_df2$Gender <- add.sample.size(plot_df2$Gender)
plot_df2$`Diagnosis Type` <- add.sample.size(plot_df2$`Diagnosis Type`)
plot_df2$Treatment <- add.sample.size(plot_df2$Treatment)
plot_df2$`Tumor Location` <- add.sample.size(plot_df2$`Tumor Location`)
plot_df2$Diagnosis <- add.sample.size(plot_df2$Diagnosis)

# 4. Define color schemes
col_Diagnosis <- c("#5fff57","mediumorchid2","darkgreen","#f268d6","lightseagreen","#005082")
names(col_Diagnosis) <- names(table(plot_df2$Diagnosis))

col_gender <- c("#CC0303","#0707CF")
names(col_gender) <- names(table(plot_df2$Gender))

col_age <- c("#00F800","#146CF6","white")
names(col_age) <- names(table(plot_df2$Age))

col_dtype <- c("#cee397","#827397","#363062","#005082")
names(col_dtype) <- names(table(plot_df2$`Diagnosis Type`))

col_annot <- c("black","lightgray")
names(col_annot) <- names(table(plot_df2$Treatment))

col_tumor_loc <- c("#D4806C","#94004C","#344C68","#7C8F97")
names(col_tumor_loc) <- names(table(plot_df2$`Tumor Location`))

col_fun1 <- list(c(col_Diagnosis, col_gender, col_age, col_dtype, col_annot, col_tumor_loc))[[1]] %>%
  unlist() %>% as.list()

# Color mapping for continuous age variable
col_fun_age <- colorRamp2(breaks = summary(plot_df$Age_at_Initial_Diagnosis),
                          colors = c("#146CF6", "#188AF0", "#00B7D8", "#00D4B0", "#00E54B", "#00F800"))

col_fun2 <- col_fun1
col_fun2$`[0]` <- 'white'

col.bg <- rep('white', 3)
names(col.bg) <- unique(plot_df2 %>% dplyr::select(Age))

# 5. Generate circular heatmap
circos.clear()
pdf(file = file.path(output_dir, "Figure1A_clinical_overview.pdf"), width = 10, height = 10)

circos.par(start.degree = 90+(360-5*3-89-k)/(92+k)*k/2+0.5+k/2-1, 
           gap.degree = 1, gap.after = c(5), points.overflow.warning = FALSE)

# Age track
circos.heatmap(mat = plot_df2 %>% dplyr::select(Age), 
               split = split2,
               col = col.bg, track.height = 0.001, 
               cell.border = "white", bg.lwd = 2, cell.lwd = 2)

# Clinical annotation tracks
circos.heatmap(plot_df2 %>% dplyr::select(c('Diagnosis Type','Tumor Location','Treatment','Diagnosis','Gender')), 
               col = unlist(col_fun2),
               track.height = 0.3, 
               cell.border = "white", bg.lwd = 2, cell.lwd = 2)

# Continuous age track
circos.heatmap(mat = plot_df2 %>% dplyr::select(Age_at_Initial_Diagnosis), 
               col = col_fun_age, track.height = 0.06, 
               cell.border = "white", bg.lwd = 2, cell.lwd = 2)

# Labels
circos.labels(rep(unique(split2)[3:2], c(4,3)), 
              x = c(0,14, 37,59, 18, 26,34), 
              labels = c('0yr','5yr','10yr','15yr','20yr','30yr','40yr'),
              connection_height = mm_h(1), cex = 2/3)

highlight.sector("[0]", col = "white")

# Axis labels
cex.name <- 1.2
facing.name <- "bending.outside"

circos.trackText(x = k/2, y = 0.5,
                 labels = 'Age', cex = cex.name, 
                 sectors = '[0]', track = 3, facing = facing.name, niceFacing = TRUE)

circos.trackText(x = rep(k/2,5), y = 4:0+0.5,
                 labels = c('Diagnosis Type','Tumor Location','Treatment','Diagnosis','Sex'),
                 cex = cex.name, sectors = rep('[0]',5), track = 2, 
                 facing = facing.name, niceFacing = TRUE)

# 6. Add legends
lgd_diagnosis <- Legend(title = "Diagnosis", at = names(col_Diagnosis), legend_gp = gpar(fill = col_Diagnosis))
lgd_gender <- Legend(title = "Sex", at = names(col_gender), legend_gp = gpar(fill = col_gender))
lgd_age <- Legend(title = "Age (years)", col_fun = col_fun_age, at = c(0, 15, 39.9), labels = c(0, 15, 39.9))
lgd_dtype <- Legend(title = "Diagnosis Type", at = names(col_dtype), legend_gp = gpar(fill = col_dtype))
lgd_annot <- Legend(title = "Treatment", at = names(col_annot), legend_gp = gpar(fill = col_annot))
lgd_tumor_location <- Legend(title = "Tumor Location", at = names(col_tumor_loc), legend_gp = gpar(fill = col_tumor_loc))

# Legend layout
gap.temp <- lgd_diagnosis@grob$vp$width - lgd_tumor_location@grob$vp$width - lgd_dtype@grob$vp$width

h0 <- lgd_age@grob$vp$height
h1 <- lgd_diagnosis@grob$vp$height
h2 <- lgd_tumor_location@grob$vp$height
h3 <- lgd_annot@grob$vp$height

h <- dev.size()[2]

circle_size <- unit(1, "snpc")

lgd_list <- packLegend(lgd_age, lgd_gender, max_height = unit(0.9*h, "inch"), direction = "horizontal",
                       column_gap = unit(0.5, "inch"))
lgd_list1 <- packLegend(lgd_diagnosis, max_height = unit(0.9*h, "inch"), direction = "horizontal")
lgd_list2 <- packLegend(lgd_tumor_location, lgd_dtype, max_height = unit(0.9*h, "inch"), direction = "horizontal",
                        column_gap = gap.temp)
lgd_list3 <- packLegend(lgd_annot, max_height = unit(0.9*h, "inch"), direction = "horizontal")

lgd.x <- unit(h/2, "inch")
lgd.gap <- unit(5, "mm")
lgd.y3 <- (unit(254, "mm")-(lgd.gap*3+h0+h1+h2))/2

draw(lgd_list,  x = lgd.x, y = lgd.y3+lgd.gap*3+h0/2+h1+h2+h3/2) 
draw(lgd_list1, x = lgd.x, y = lgd.y3+lgd.gap*2+h1/2+h2+h3/2) 
draw(lgd_list2, x = lgd.x, y = lgd.y3+lgd.gap*1+h2/2+h3/2) 
draw(lgd_list3, x = lgd.x, y = lgd.y3) 

# 7. Add text labels for cohorts
cex.n <- 0.8
x.n <- -0.1

text(x = x.n, y = 0.38, 'AYA', cex = cex.n)
text(x = x.n, y = 0.34, '(n = 34)', cex = cex.n)

text(x = x.n, y = 0.28, 'PED', cex = cex.n)
text(x = x.n, y = 0.24, '(n = 59)', cex = cex.n)

# 8. Close plot
circos.clear()
dev.off()
