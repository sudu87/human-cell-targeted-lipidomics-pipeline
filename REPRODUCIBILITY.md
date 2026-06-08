# Reproducibility Instructions

This document gives a minimal, repeatable way to rerun the analyses in this repository.

## 1. Clone The Repository

```sh
git clone https://github.com/sudu87/human-cell-targeted-lipidomics-pipeline.git
cd human-cell-targeted-lipidomics-pipeline
```

## 2. Restore The R Environment

The repository includes `renv.lock`, which records the R and package versions used for the scripts.

```r
install.packages("renv")
renv::restore()
```

The lockfile was generated with R 4.5.2. If you use another R version, `renv` may still restore the package environment, but exact binary package builds can differ by operating system and R version.

## 3. Download Input Data

Input data are deposited in Zenodo:

Das, S., & Mohanty, A. (2026). Targeted lipidomics analysis of sphingolipid metabolism during Simkania negevensis infection (v.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.18866967

Download the Zenodo files into a local `data/` directory:

```text
human-cell-targeted-lipidomics-pipeline/
├─ data/
├─ outputs/
├─ scripts/
├─ documents and instructions/
├─ README.md
└─ renv.lock
```

The `data/` and `outputs/` directories are local working directories. Input data are not tracked in this GitHub repository.

## 4. Configure Script Paths

Each R script is standalone. Before running a script, open it and edit the user-defined paths near the top.

Common fields are:

```r
input_file <- "data/your_input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "outputs/your_analysis_name"
```

Some scripts use `file_path`, `out_dir`, or an `analysis_configs` list instead of `input_file` and `output_dir`. Edit those fields in the same way.

Use exact sheet names from the Excel workbooks. The scripts apply `janitor::clean_names()` to column names, but they do not change sheet names.

## 5. Suggested Data-To-Script Mapping

These are the main Zenodo files expected by the current scripts.

| Analysis | Zenodo input file | Main script(s) |
| --- | --- | --- |
| HeLa/HUVEC total sphingolipids | `1. Raw data_HUVEC_Hela-Total Sphingolipids.xlsx` | `pca_permanova_hela_huvec_analysis.R`, `cell_line_averaged_lipid_heatmap_analysis.R` |
| HeLa/HUVEC individual sphingolipids | `2. Raw data_HUVEC_Hela-individual sphingolipids.xlsx` | `cell_line_averaged_lipid_heatmap_analysis.R` |
| SMS1/2 total lipids | `29. Raw data_R_SMS1&2_TotalSL.xlsx` | `sms12_manova_pca_dispersion_analysis.R`, `sms12_lipid_summary_barplots.R`, `sms12_lipid_composition_stacked_barplot.R`, `sms12_ordered_condition_heatmap.R`, `sms12_per_lipid_anova_posthoc_analysis.R`, `transformation_diagnostics.R` |
| SMS1/2 individual lipids | `31. Raw data_R_SMS1&2_IndividualSL.xlsx` | `sms12_individual_lipids_permanova_dispersion_analysis.R`, `sms12_per_lipid_anova_posthoc_analysis.R`, `transformation_diagnostics.R` |
| SMS1/2 total pairwise contrasts | `30. Pairwise comparasion_ SMS1:2_ Total SL.xlsx` | `custom_pairwise_lipid_contrast_heatmap.R` |
| SMS1/2 individual pairwise contrasts | `32. Pairwise comparasion_ SMS1:2_ IndividualSL.xlsx` | `custom_pairwise_lipid_contrast_heatmap.R` |
| Inhibitor total lipid raw data | `5. Rawdata_ARC39_Total SL.xlsx`, `6. Rawdata_Desipramin_Total SL.xlsx`, `7. Rawdata_HPA-12_Total SL.xlsx`, `8. Rawdata_Myriocin_Total SL.xlsx`, `9. Rawdata_AKS466_Total SL.xlsx`, `10. Rawdata_untreated_Total SL.xlsx`, `11. Rawdata_DMSO_Total SL.xlsx` | `pca_permanova_inhibitor_sphingolipids_analysis.R`, `inhibitor_infection_heatmap_average_replicates.R`, `manova_per_lipid_anova_analysis.R` |
| Inhibitor individual lipid raw data | `17. Rawdata_ARC39_Individual SL.xlsx`, `18. Rawdata_AKS466_Individual SL.xlsx`, `19. Rawdata_Desipramine_Individual SL.xlsx`, `20. Rawdata_HPA-12_Individual SL.xlsx`, `21. Rawdata_Myriocin_Individual SL_new.xlsx`, `22. Rawdata_Untreated_Individual SL.xlsx`, `23. Rawdata_DMSO_Individual SL.xlsx` | `inhibitor_infection_heatmap_average_replicates.R`, `manova_per_lipid_anova_analysis.R` |
| Inhibitor total pairwise contrasts | `12. Pairwise comparasion_ARC39_posthoc inhibitor contrasts by infection_Total SL.xlsx`, `13. Pairwise comparasion_Desipramine_posthoc inhibitor contrasts by infection_Total SL.xlsx`, `14. Pairwise comparasion_HPA-12_posthoc inhibitor contrasts by infection_Total SL.xlsx`, `15. Pairwise comparasion_Myriocin_posthoc inhibitor contrasts by infection_Total SL.xlsx`, `16. Pairwise comparasion_AKS466_posthoc inhibitor contrasts by infection_Total SL.xlsx` | `pairwise_fdr_lipid_heatmap.R`, `custom_pairwise_lipid_contrast_heatmap.R` |
| Inhibitor individual pairwise contrasts | `24. Pairwise comparasion_ARC39_posthoc inhibitor contrasts by infection_Individual SL.xlsx`, `25. Pairwise comparasion_Desipramine_posthoc inhibitor contrasts by infection_Individual SL.xlsx`, `26. Pairwise comparasion_HPA-12_posthoc inhibitor contrasts by infection_Individual SL.xlsx`, `27. Pairwise comparasion_Myriocin_posthoc inhibitor contrasts by infection_Individual SL.xlsx`, `28. Pairwise comparasion_AKS466_posthoc inhibitor contrasts by infection_Individual SL.xlsx` | `pairwise_fdr_lipid_heatmap.R`, `custom_pairwise_lipid_contrast_heatmap.R` |

## 6. Run Scripts

Run scripts from the repository root. For example:

```sh
Rscript scripts/transformation_diagnostics.R
Rscript scripts/sms12_manova_pca_dispersion_analysis.R
Rscript scripts/sms12_lipid_summary_barplots.R
Rscript scripts/sms12_lipid_composition_stacked_barplot.R
Rscript scripts/sms12_ordered_condition_heatmap.R
Rscript scripts/sms12_individual_lipids_permanova_dispersion_analysis.R
Rscript scripts/custom_pairwise_lipid_contrast_heatmap.R
```

There is no single mandatory order for every script. Most scripts start from raw Excel data or pairwise contrast tables and can be run independently once their input path, sheet name, and output directory are set.

## 7. Recommended SMS1/2 Run

For the SMS1/2 total lipid workflow, configure these scripts with:

```r
input_file <- "data/29. Raw data_R_SMS1&2_TotalSL.xlsx"
output_dir <- "outputs/sms12_total"
```

Set `sheet_name` to the exact sheet containing the SMS1/2 total lipid table.

Then run:

```sh
Rscript scripts/transformation_diagnostics.R
Rscript scripts/sms12_manova_pca_dispersion_analysis.R
Rscript scripts/sms12_lipid_summary_barplots.R
Rscript scripts/sms12_lipid_composition_stacked_barplot.R
Rscript scripts/sms12_ordered_condition_heatmap.R
Rscript scripts/sms12_per_lipid_anova_posthoc_analysis.R
```

For SMS1/2 individual lipid PERMANOVA, configure:

```r
input_file <- "data/31. Raw data_R_SMS1&2_IndividualSL.xlsx"
output_dir <- "outputs/sms12_individual"
```

Set `sheet_name` to the exact sheet containing the SMS1/2 individual lipid table, then run:

```sh
Rscript scripts/sms12_individual_lipids_permanova_dispersion_analysis.R
```

For SMS1/2 pairwise contrast heatmaps, configure `custom_pairwise_lipid_contrast_heatmap.R` with either:

```r
input_file <- "data/30. Pairwise comparasion_ SMS1:2_ Total SL.xlsx"
output_dir <- "outputs/sms12_pairwise_total"
```

or:

```r
input_file <- "data/32. Pairwise comparasion_ SMS1:2_ IndividualSL.xlsx"
output_dir <- "outputs/sms12_pairwise_individual"
```

Then run:

```sh
Rscript scripts/custom_pairwise_lipid_contrast_heatmap.R
```

## 8. Check Outputs

Each script writes files to the configured output directory. Typical outputs include:

- diagnostic plots for raw and log-transformed values
- heatmaps
- PCA plots and PCA variance tables
- MANOVA or PERMANOVA text summaries
- dispersion test summaries
- per-lipid ANOVA and posthoc tables
- pairwise contrast heatmap plot data

Keep output directories separate by analysis so that files from different workflows do not overwrite each other.

## 9. Record Provenance

For a fully reproducible rerun, record:

- Git commit used for the analysis
- Zenodo DOI and file version
- R version
- `renv.lock` version
- script name
- input file name
- Excel sheet name
- output directory
- any edited condition, infection, contrast, or lipid-order vectors

The current Git commit can be recorded with:

```sh
git rev-parse HEAD
```

