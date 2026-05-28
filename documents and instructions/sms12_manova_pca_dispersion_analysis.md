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

## User-defined paths

Open `scripts/sms12_manova_pca_dispersion_analysis.R` and set the input file, sheet name, and output directory:

```r
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"
```

Use the exact sheet name from your Excel workbook for `sheet_name`.

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

Outputs are written to the user-defined `output_dir` and its `statistics` subfolder:

```r
plot_dir <- output_dir
out_dir <- file.path(output_dir, "statistics")
```

The script writes:

```text
<output_dir>/pca_variance_explained.csv
<output_dir>/pca_scores.csv
<output_dir>/pca_scree.pdf
<output_dir>/pca_scores.pdf
<output_dir>/statistics/manova_pillai.txt
<output_dir>/statistics/dispersion_test.txt
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
