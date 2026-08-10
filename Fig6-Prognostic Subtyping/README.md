# Figure 6 — Prognostic Subtyping

Complete the shared setup steps in the root `README.md` before running these scripts. Unless otherwise noted, run scripts from within this folder:

```bash
cd "Fig6-Prognostic Subtyping"
```

Generated PDFs, PNG previews, and tables are written to `output/`.

## Figure Panels

| Analysis/Figures/Tables | Script | Section in the manuscript | Simplified legend |
| :---: | :--- | :--- | :--- |
| Figure 6A | Not yet represented by a standalone public script | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Heatmap showing distinct protein pathway activities across protein clusters identified in the Discovery cohort (ages 0-62). |
| Figure 6B | Not yet represented by a standalone public script | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Heatmap illustrating protein pathway ssGSEA scores for the same pathways as shown in A within the validation data. |
| Figure 6C | Not yet represented by a standalone public script | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Bubble heatmap illustrating correlations between average pathway scores of each protein cluster in the cDiscovery cohort versus the Validation cohort. |
| Figure 6D | `Figure6D_S6B_cluster_annotation_corrplots.R` | Protein Pathway Activities Reveal Distinct Subtypes of High Grade Glioma | Bubble heatmap illustrating distributions of H3-3A mutated tumors across protein clusters. |
| Figure 6E | `Figure6M_S6O_S6E_S6J_6I_cd276_itgav_pla2g4a_subtype_analysis.R` | Causal kinases and CAR-T markers associated with Protein Subtypes | Boxplots of ATM, LCK, NTRK2, and CDK8 kinase activities in the cDiscovery cohort (ages 0-40), stratified by protein clusters. |
| Figure 6F | `Figure6M_S6O_S6E_S6J_6I_cd276_itgav_pla2g4a_subtype_analysis.R` | Causal kinases and CAR-T markers associated with Protein Subtypes | Boxplot of CD276 (B7-H3) protein abundance in the cDiscovery cohort (ages 0-40), stratified by sex and protein clusters. |
| Figure 6G | `Figure6G_subtype_survival_analysis.R` | Overall survival outcome differences across Protein Subtypes | Survival curves comparing F2, M2, and other clusters within the cDiscovery cohort (ages 0-40). |
| Figure 6H | `Figure6M_S6O_S6E_S6J_6I_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Volcano plot comparing protein abundances between F2 and M2 for genes in the Interferon Gamma Response pathway. |
| Figure 6I | `Figure6M_S6O_S6E_S6J_6I_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Boxplots illustrating protein and phosphosite abundance of PLA2G4A across sexes and protein clusters. |
| Figure 6J | `Figure6M_S6O_S6E_S6J_6I_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Scatterplot showing correlations between LCK kinase activity and PLA2G4A protein expression within each sex in the cDiscovery cohort (ages 0-40). |
| Figure 6K | `Figure6K_sex_biased_glycopeptide_pathway_enrichment.R` | Overall survival outcome differences across Protein Subtypes | Bar plot illustrating pathway enrichment results for comparing male and female tumors in C2 based on glycopeptide abundance data. |
| Figure 6L | `Figure6K_sex_biased_glycopeptide_pathway_enrichment.R` | Overall survival outcome differences across Protein Subtypes | Volcano plot comparing glycopeptide abundances between F2 and M2 for genes in the Epithelial Mesenchymal Transition pathway. |
| Figure 6M | `Figure6M_S6O_S6E_S6J_6I_cd276_itgav_pla2g4a_subtype_analysis.R` | Overall survival outcome differences across Protein Subtypes | Boxplot illustrating abundance of the ITGAV glycopeptide with an N-linked sialylated glycan (N3H5F1S1G0) across sex and protein clusters. |
| Figure 6N | `Figure6N_itgav_sialylated_glycopeptide_pdl1_correlation.R` | Overall survival outcome differences across Protein Subtypes | Scatterplot showing correlations between global protein abundance of PD-L1 and abundance of the ITGAV glycopeptide N3H5F1S1G0 in female and male tumors. |

Note: `Figure6D_S6B_cluster_annotation_corrplots.R` generates Figure 6D and
Figure S6B annotation-enrichment corrplots. It does not generate the Figure 6A,
Figure 6B, or Figure 6C protein pathway activity heatmaps/correlation panel.
