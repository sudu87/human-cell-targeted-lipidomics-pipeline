# PCA and PERMANOVA analysis for total sphingolipids

This document explains how to use `scripts/pca_permanova_total_sphingolipids_analysis.R`.

The pipeline is designed for total sphingolipid Excel files where each file represents one inhibitor or control condition. It keeps the two solvent/control groups separate:

- DMSO-soluble inhibitor group, analyzed with the DMSO control
- Water-soluble inhibitor group, analyzed with the untreated control

Do not combine the DMSO and untreated groups into one analysis unless you explicitly model solvent/control group as a separate variable.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
stringr
purrr
ggplot2
tibble
vegan
tools
```

Install missing packages before running the script.

## Input folder structure

The script expects the following input folders by default:

```text
averaged_totalsphingolipids/new files/with_DMSO/
averaged_totalsphingolipids/new files/with untreated/
```

Each folder should contain `.xlsx` files. Each Excel file should correspond to one inhibitor or control condition.

For the DMSO-soluble group, the default file-name mapping is:

```r
name_map = c(
  "AKS466_only_Total SL" = "AKS466 total",
  "DMSO_only_Total SL" = "DMSO total",
  "HPA-12_only_Total SL" = "HPA-12 total",
  "Myriocin_only_Total SL" = "Myriocin total"
)
```

For the water-soluble/untreated group, the default file-name mapping is:

```r
name_map = c(
  "ARC39_only_Total SL" = "ARC39 total",
  "untreated_only_Total SL" = "untreated total",
  "Desipramin_only_Total SL" = "Desipramin total"
)
```

The names on the left must match the Excel filenames without `.xlsx`.

## Input data requirements

Each Excel file must:

- contain one row per sample
- contain a metadata column named `infection`
- contain lipid measurements in the remaining columns
- use numeric lipid values, or values that can be converted to numeric
- represent one inhibitor/control condition per workbook

The script adds the `inhibitor` column automatically based on the filename.

## Choosing which analysis to run

Open `scripts/pca_permanova_total_sphingolipids_analysis.R` and set:

```r
ANALYSIS_GROUP <- "dmso"
```

to analyze the DMSO-soluble group.

Use:

```r
ANALYSIS_GROUP <- "untreated"
```

to analyze the water-soluble group with untreated control.

## What the script does

The script performs the following steps:

- reads all `.xlsx` files from the selected input folder
- assigns inhibitor/control names using the selected `name_map`
- stops with an error if an unexpected filename is found
- cleans column names with `janitor::clean_names()`
- aligns files that have different lipid columns
- combines all files into one data frame
- converts `infection` and `inhibitor` into factors
- identifies lipid columns as all columns except `infection` and `inhibitor`
- coerces lipid columns to numeric
- applies `log10(x + 1)` transformation
- removes lipid columns with only missing values
- removes lipid columns with zero variance
- z-scores lipid columns when `METHOD <- "euclidean"`
- removes incomplete sample rows
- runs PCA on the filtered lipid matrix
- saves PCA variance explained
- creates and saves PCA plots
- runs PERMANOVA testing the effect of `inhibitor`
- runs dispersion checks using `betadisper()` and `permutest()`

## Output files

For the DMSO group, outputs are written by default to:

```text
averaged_totalsphingolipids/new files/with_DMSO/with_dmso_analysis_outputs_Mar2026/
```

For the untreated group, outputs are written by default to:

```text
averaged_totalsphingolipids/new files/with untreated/with_untreated_analysis_outputs_Mar2026/
```

The script writes:

```text
<analysis_group>_pca_variance_explained.csv
<analysis_group>_pca_plot.pdf
<analysis_group>_pca_plot.png
<analysis_group>_permanova_main.txt
<analysis_group>_permanova_dispersion_tests.txt
```

For example, when `ANALYSIS_GROUP <- "dmso"`, the files are:

```text
dmso_pca_variance_explained.csv
dmso_pca_plot.pdf
dmso_pca_plot.png
dmso_permanova_main.txt
dmso_permanova_dispersion_tests.txt
```

## Main R objects created

The script creates:

- `df`: combined data before complete-case filtering
- `df_cc`: complete-case metadata and transformed lipid data
- `X_cc`: complete-case processed lipid matrix
- `pca_fit`: PCA model from `prcomp()`
- `pca_var_df`: variance explained by each principal component
- `scores`: PCA scores with metadata
- `p_pca`: PCA plot object
- `adon_inhibitor`: PERMANOVA result
- `per_inhib`: dispersion test result

## Notes

- `METHOD <- "euclidean"` is appropriate because the lipid matrix is z-scored before PCA and PERMANOVA.
- The script uses complete-case filtering after QC and scaling, so samples with missing lipid values are removed.
- `MIN_N <- 4` prevents analysis when too few complete observations remain.
- If a filename is not present in `name_map`, the script stops instead of assigning `NA` as the inhibitor.
- The DMSO group and untreated group are analyzed separately because their controls are biologically different.
