# Figure 4 — Prognostic Markers and Pathways Based on AD-TMP

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig4-Prognostic markers:pathways based on AD-TMP"
```

Generated PDFs, PNG previews, model objects, and tables are written to `output/`.

Scripts that call `library(temporalCPSA)` or `temporalCPSA::` require installation of `temporalCPSA` from GitHub as described in the root `README.md`. In this folder, `Figure4_sex_specific_cox_models.R` is the only script that directly requires `temporalCPSA`; `Figure4_combined_frozen_score_models.R` uses the `Figure4_analysis_ready_data.rds` file produced by that script. If the package source has changed, reinstall `temporalCPSA` before running `Figure4_sex_specific_cox_models.R`.

## Figure Panels

| Analysis/Figures/Tables | Script | Main input files | Primary output files |
| :--- | :--- | :--- | :--- |
| Figure 4A | `Figure4_sex_specific_cox_models.R` | `../data/STable1.xlsx`, `../data/STable4.xlsx`, `../data/cDisc_mutation_10192023.tsv`, `../data/cDisc_proteome_imputed_data_09152023.tsv` | This analysis is part of the sex-specific prognostic modeling workflow. The current script emits the named Figure 4/S4 model outputs listed below rather than a separate `Figure4A_*.pdf` file. |
| Figure 4B | `Figure4B_survival_assoc_glyco.R` | `../data/STable4.xlsx`, sheets `SA-Glyco-Disc` and `SA-Glyco-Pathway-cDisc` | `Figure4B_survival_assoc_glyco.pdf`; `Figure4B_survival_assoc_glyco.png`; `Figure4B_selected_glycopeptides_from_STable4.tsv` |
| Figure 4C | `Figure4C_survival_landscape.R` | `../data/STable4.xlsx`, sheets `SA-Protein-cDisc-Ref` and `SA-RNA-cDisc-Ref` | `Figure4C_survival_landscape.pdf` |
| Figure 4D | `Figure4D_survival_assoc_heatmap.R` | `../data/STable4.xlsx`, sheets `SA-Protein-cDisc-Ref` and `SA-Protein-Pathway-cDisc-Ref` | `Figure4D_protein_survival_overlap.pdf` |
| Figure 4E; Figure S4D; Figure S4E; Figure S4F; Figure 4I; Figure S4I; STable4 model tables | `Figure4_sex_specific_cox_models.R` | `../data/STable1.xlsx`, `../data/STable4.xlsx`, `../data/cDisc_mutation_10192023.tsv`, `../data/cDisc_proteome_imputed_data_09152023.tsv` | `Figure4E_weights_sel_cDiscovery_Male.pdf`; `FigureS4D_prog_combined_Male_cDiscovery.pdf`; `FigureS4E_prog_combined_Female_cDiscovery.pdf`; `FigureS4F_weights_sel_cDiscovery_Female.pdf`; `Figure4I_PHKG2_Male_cDiscovery_splot.pdf`; `FigureS4I_MORC2_Male_cDiscovery_splot.pdf`; `STable4_multivariate_hr_table_Male_cDiscovery.tsv`; `STable4_multivariate_hr_table_Female_cDiscovery.tsv`; `Figure4_analysis_ready_data.rds` |
| Figure 4G; Figure 4H; Figure S4 combined-score panels | `Figure4_combined_frozen_score_models.R` | `Figure4_analysis_ready_data.rds` produced by `Figure4_sex_specific_cox_models.R` | `Figure4G_combined_score_HR_Both_cDiscovery.pdf`; `Figure4H_combined_score_HR_Both_Validation.pdf`; `FigureS4F_combined_male_score_cDiscovery.pdf`; `FigureS4F_combined_female_score_cDiscovery.pdf`; `FigureS4_combined_male_score_Validation.pdf`; `FigureS4_combined_female_score_Validation.pdf`; `multivariate_hr_table_Both_cDiscovery.tsv`; `multivariate_hr_table_Both_Validation.tsv` |
| Figure 4J and related validation outputs | `Figure4_validation_combined_score_model.R` | `../data/STable1.xlsx`, validation clinical and protein sheets | `MORC2_MaleDerived_Validation_hr_plot.pdf`; `PHKG2_MaleDerived_Validation_hr_plot.pdf`; `MORC2_MaleDerived_Validation_splot.pdf`; `PHKG2_MaleDerived_Validation_splot.pdf`; `MORC2_MaleDerived_Validation_hr_table.tsv`; `PHKG2_MaleDerived_Validation_hr_table.tsv` |

## Run Order

The table-only panels can be run independently:

```bash
Rscript Figure4B_survival_assoc_glyco.R
Rscript Figure4C_survival_landscape.R
Rscript Figure4D_survival_assoc_heatmap.R
```

For Figure 4B, pathway groups are taken from the five rows marked
`selected == "Yes"` in `SA-Glyco-Pathway-cDisc`. Glycopeptide rows are then
selected from those pathway/module gene sets using the original glyco
`glyco.int2` PED FDR criterion and the strongest significant membrane
glycopeptide per gene within sex before union. Glycopeptide-level heatmap
values and FDR point overlays are read from `SA-Glyco-Disc`. `FAM234B` passes
this reconstructed `CLUSTER_24` rule but was not present in the original
displayed Figure 4B panel; Table 4 does not encode a final displayed-panel flag,
so the script documents this discrepancy rather than treating `FAM234B` as a
failed candidate.

The prognostic-score analyses should be run in this order:

```bash
Rscript Figure4_sex_specific_cox_models.R
Rscript Figure4_combined_frozen_score_models.R
Rscript Figure4_validation_combined_score_model.R
```

`Figure4_sex_specific_cox_models.R` writes `Figure4_analysis_ready_data.rds`, which stores the harmonized cDiscovery and Validation objects plus frozen score coefficients used by `Figure4_combined_frozen_score_models.R`.
