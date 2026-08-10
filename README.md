# Pediatric/AYA High-Grade Glioma (HOPE Study) — Analysis Repository

This repository contains code used to generate all figures and tables for the HOPE pediatric and adolescent/young adult (AYA) high-grade glioma study.

---

# Reproducibility Overview

1. Download supplementary tables from the journal website  
2. Download raw data from the PDC  
3. (Optional) Install `temporalCPSA` for reusable trajectory/CPSA analyses  
4. Follow the README for the figure of interest  

---

# Repository Structure

Each folder corresponds to a figure in the manuscript and contains the scripts required to reproduce that analysis.

## Figure-to-Folder Mapping

| Manuscript Figure | Repository Folder |
|---|---|
| Figure 1 | `Fig1-Study overview and Landscape of DNA alteration` |
| Figure 2 | `Fig2-Age-dependent tumor molecular profiles` |
| Figure 3 | `Fig3-Cis and Trans Regulation` |
| Figure 4 | `Fig4-Prognostic markers:pathways based on AD-TMP` |
| Figure 5 | `Fig5-Causal Kinase Analysis` |
| Figure 6 | `Fig6-Prognostic Subtyping` |
| Figure 7 | `Fig7-Immune Landscape` |

---

# Data Availability and Setup

## Step 1 — Download Supplementary Tables

Download Supplementary Tables `STable1`–`STable7` associated with the manuscript from the journal website.

Place all seven supplementary table files directly into:

```text
data/
```

Important:

- Do not create subdirectories
- Do not rename files

---

## Step 2 — Download Raw Data from the PDC

1. Go to the PDC website:

   ```text
   https://proteomic.datacommons.cancer.gov/pdc/
   ```

2. Search for one of the following study IDs:

   ```text
   PDC000497
   PDC000498
   PDC000499
   ```

3. Open the study summary page.

4. In the **Supplementary Data** panel:

   - Locate the row labeled:

     ```text
     Publication Supplementary Material (Archive)
     ```

   - Click the file count, 1

5. Download:

   ```text
   HOPE_AYA_supplementary_files_updated_05032026.zip
   ```

6. Unzip the archive.

The archive contains the 14 processed data files required for analyses associated with:

```text
PDC000497
PDC000498
PDC000499
```

Place all extracted files directly into:

```text
data/
```

Important:

- Do not create subdirectories
- Do not rename files

---

# Step 3 (Optional) — Install temporalCPSA

Reusable trajectory and Cross-Population Survival Analysis (CPSA) methods are
maintained in a separate R package:

[`temporalCPSA`](https://github.com/lashbroz/temporalCPSA)

Package documentation and the overview article are available at:
<https://lashbroz.github.io/temporalCPSA/>

`temporalCPSA` provides the reusable analysis framework for age-dependent tumor
molecular trajectory analyses, tumor-normal/reference trajectory comparisons,
trajectory clustering, optional exploratory age-class diagnostics, and
downstream CPSA. The package is intended for both reproducible manuscript
analyses and reuse in independent cohorts.

Install `temporalCPSA` from GitHub before running updated figure scripts that
call package functions:

```r
install.packages("remotes")
remotes::install_github("lashbroz/temporalCPSA")
```

Then confirm installation:

```r
library(temporalCPSA)
temporalCPSA::ageTMP_status()
```

Important:

- `temporalCPSA` is required only for scripts that call package functions such
  as `ageTMP_predict_tumor_trajectory_matrix()` or
  `ageTMP_fit_reference_cpsa()`.
- The package is designed to read from documented files in `data/` and from
  documented package reference data where applicable.

---

# Step 4 — Run Figure-Specific Analyses

1. Navigate to the folder corresponding to the figure of interest  
2. Open that folder’s `README.md`  
3. Run the scripts described there to reproduce the analysis  

Each figure folder contains self-contained code and any figure-specific instructions.

When scripts are run from within a figure folder, shared manuscript data are
read from:

```text
../data/
```

Generated figure files and tables are written to the `output/` subdirectory
within the respective figure folder, unless a figure-specific README explicitly
documents an exception.

---

# Requirements

- Analyses were developed in R
- Required R packages are loaded within each script
- Scripts assume all required input files are present in the `data/` directory
