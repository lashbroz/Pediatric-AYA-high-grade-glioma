# Figure 1 — Study Overview and Landscape of DNA Alteration

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig1-Study overview and Landscape of DNA alteration"
```

Generated PDFs, PNG previews, and tables are written to `output/`.

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Figure 1A | `Figure1A_clinical_overview.R` | Proteogenomics profiling of the Discovery cohort | Summary of the demographic and clinical annotations of the Discovery study cohort (N=93). |
| Figure 1B | `Figure1B_barplot.R` | Methylation and Genomically informed Cancer Subtype Groups | Distribution of cancer groups stratified by PED and AYA. |
| Figure 1C | `figure1c_oncplot.R` | Age dependent Mutational and CNV Landscape of HGG | Oncoprint plot for the frequent (>5%) mutations detected in the Discovery Cohort (N=89). |
| Figure 1D | `Figure1E_reference_assoc.R` | Age dependent Mutational and CNV Landscape of HGG | Mutation frequency profiles along age for two groups of genes with distinct mutation patterns based on the Reference cohort data. |
| Figure 1E | `Figure1E_reference_assoc.R` | Age dependent Mutational and CNV Landscape of HGG | Association of overall survival with somatic mutation events across age groups in the Reference cohort, evaluated using multivariate Cox regression models. |
| Figure 1F | `Figure1F_AD_TMP_heatmap.R`; `Figure1F_AD_TMP_heatmap_source_derived.R`; `Figure1F_AD_TMP_age_segmentation_support.R` | Transcriptomic and Proteomic Profiles Defining Adolescent versus Young Adult Age Groups | The first script reproduces the manuscript panel using the archived manuscript visualization input; the second renders a source-derived analogue from the final repository data tables using `temporalCPSA`; the third summarizes data-driven contiguous age-segmentation support from the source-derived AD-TMP matrix. |
| Figure 1G | `Figure1G_reference_survival_assoc.R` | Transcriptomic and Proteomic Profiles Defining Adolescent versus Young Adult Age Groups | Kaplan-Meier overall survival curves for different age groups in the Reference cohort. |
| Figure S1A | `FigureS1A_clinical_overview.R` | Supplementary cohort annotations | Related cohort annotation summary. |
| Figure S1D | `FigureS1D_CNV_landscape.R` | Age dependent Mutational and CNV Landscape of HGG | Copy-number landscape overview for the Discovery cohort. |

## Run Order

Most Figure 1 scripts can be run independently.

### Figure 1F

Figure 1F has companion scripts for manuscript visualization reproduction,
source-derived reanalysis, and data-driven age-segmentation support. The
manuscript reproduction script renders the published panel using the archived
visualization input generated during figure preparation. The Figure 1F scripts
require `temporalCPSA`.

```bash
Rscript Figure1F_AD_TMP_heatmap.R
```

The age-class structure shown in Figure 1F was derived early in the study using
the data versions available at that stage of manuscript development. To show the
same visualization concept on the final repository data tables, run the
source-derived analogue:

```bash
Rscript Figure1F_AD_TMP_heatmap_source_derived.R
```

The source-derived analogue uses `temporalCPSA` functions to generate the
AD-TMP matrix from the final source tables before plotting. It is intended to
illustrate the age-class distinctions in the final dataset and is not expected
to be pixel-identical to the manuscript panel. The main differences are data
version and feature-universe differences rather than a different sample-ordering
scheme. Protein trajectories are reproduced with near-exact agreement. RNA uses
the final transcriptomic source matrix, whose gene and sample universe differs
modestly from the earlier matrix used during figure development. Phospho is
estimated at the phosphosite level and then collapsed to `ApprovedGeneSymbol`;
the final source phospho table contains a substantially larger phosphosite
universe, so a given gene can be represented by a different set of phosphosites
than in the archived visualization input.

The survival-days line in Figure 1F is a descriptive annotation rather than a
formal survival model. The manuscript reproduction uses the archived
figure-preparation survival annotation, whereas the source-derived analogue
defaults to the currently available reference clinical table. The source-derived
reference subset has a lower observed-survival distribution among event-coded
cases, so its LOESS peak is lower than in the archived visualization. We retain
the same smoothing approach rather than tuning the span to visually match the
published panel.

After generating the source-derived AD-TMP matrix, run the contiguous
age-segmentation diagnostic:

```bash
Rscript Figure1F_AD_TMP_age_segmentation_support.R
```

This support script applies the exploratory `temporalCPSA` age-contiguous
segmentation diagnostic to the source-derived Figure 1F AD-TMP matrix, writes
candidate cutpoints and mean silhouette summaries, and renders the green
diagnostic plot used to visualize support for molecularly coherent age
intervals.

These scripts read manuscript data from `../data/` and write generated files to
`output/`.
