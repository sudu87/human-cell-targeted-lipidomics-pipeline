# SMS1/2 total lipid MANOVA, PCA, and dispersion analysis

This document explains how to use `scripts/sms12_manova_pca_dispersion_analysis.R`.

The script performs the first part of the SMS1/2 total lipid analysis: data import, log transformation, MANOVA, PCA, and dispersion checks.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
ggplot2
tibble
vegan
```

Install missing packages before running the script.

## Input file

The script expects the following Excel file by default:

```text
sms12_lipidomics/SMS1_2_Total lipids_no_19_test_file.xlsx
```

It reads the sheet:

```text
SMS1&2_Total lipids_uninfected
```

## Input data requirements

The Excel sheet must contain:

- one row per sample
- a metadata column named `infection`
- a metadata column named `condition`
- lipid abundance measurements in the remaining columns

The `infection` column is converted to a factor with the levels:

```r
c("no", "yes")
```

The `condition` column is converted to a factor using the condition labels present in the data.

## What the script does

The script performs the following steps:

- reads the Excel sheet with `readxl::read_excel()`
- cleans column names with `janitor::clean_names()`
- identifies lipid columns as all columns except `infection` and `condition`
- coerces lipid columns to numeric
- applies `log10(x + 1)` transformation
- keeps complete cases for the lipid matrix and metadata
- runs MANOVA using `infection * condition`
- saves Pillai's trace MANOVA output
- runs PCA using centered and scaled lipid values
- saves PCA variance explained
- saves PCA scores
- creates a PCA scree plot
- creates a PCA score plot colored by `condition` and shaped by `infection`
- tests group dispersion with `vegan::betadisper()` and `vegan::permutest()`

## Output files

Outputs are written to:

```text
sms12_lipidomics/
sms12_lipidomics/statistics/
```

The script writes:

```text
sms12_lipidomics/PCA_variance_explained_SMS1_2_total_lipids_test_April2026.csv
sms12_lipidomics/PCA_scores_SMS1_2_total_lipids_test_April2026.csv
sms12_lipidomics/pca_scree_SMS1_2_total_lipids_test_April2026.pdf
sms12_lipidomics/pca_scores_SMS1_2_total_lipids_test_April2026.pdf
sms12_lipidomics/statistics/manova_pillai_SMS1_2_total_lipids_test_April2026.txt
sms12_lipidomics/statistics/dispersion_test_SMS1_2_total_lipids_test_April2026.txt
```

## Main R objects created

The script creates:

- `sms_test`: cleaned and log-transformed input data
- `lipid_cols`: lipid measurement columns
- `resp_mat_cc`: complete-case lipid matrix
- `sms_test_cc`: complete-case metadata and transformed lipid data
- `sms_test_manova`: MANOVA model
- `manova_pillai`: Pillai's trace MANOVA summary
- `pca_sms_test`: PCA model from `prcomp()`
- `pca_var_test_df`: PCA variance explained table
- `scores_sms_test`: PCA scores with metadata
- `p_scree_sms12_test`: PCA scree plot
- `p_pca_sms12_test`: PCA score plot
- `disp_groups_test`: betadisper model
- `per_groups_test`: permutation test for dispersion

## Notes

- The script uses complete-case filtering before MANOVA, PCA, and dispersion testing.
- PCA is performed with `center = TRUE` and `scale. = TRUE`, so each lipid contributes on a comparable scale.
- Dispersion is tested on the log-transformed lipid matrix using Euclidean distance.
- The scree plot save call uses `p_scree_sms12_test`, matching the object created in the script.
