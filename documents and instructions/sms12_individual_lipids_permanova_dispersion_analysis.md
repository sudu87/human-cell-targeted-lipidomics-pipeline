# SMS1/2 individual lipid PERMANOVA and dispersion analysis

This document explains how to use `scripts/sms12_individual_lipids_permanova_dispersion_analysis.R`.

The script performs multivariate testing for SMS1/2 individual lipid species. It uses PERMANOVA rather than classical MANOVA because individual-lipid datasets often contain many lipid variables relative to the number of samples. In that setting, the MANOVA response matrix can become rank-deficient or otherwise fail matrix constraints required for classical MANOVA. PERMANOVA is used here as a distance-based multivariate alternative.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
vegan
```

Install missing packages before running the script.

## User-defined paths

Open `scripts/sms12_individual_lipids_permanova_dispersion_analysis.R` and set:

```r
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"
```

The input file should be the raw Excel workbook containing individual lipid measurements.

## Input data requirements

After `janitor::clean_names()` is applied, the Excel sheet must contain:

- one row per sample
- a column named `infection`
- a column named `condition`
- individual lipid measurement columns in the remaining columns

The script defines metadata columns as:

```r
metadata_cols <- c("infection", "condition")
```

All other columns are treated as lipid measurements. If your file contains other metadata columns, add them to `metadata_cols` before running the script.

## Analysis options

The main configurable options are:

```r
infection_order <- c("no", "yes")
method <- "euclidean"
permutations <- 999
set_seed <- 1
log_transform <- TRUE
```

The default transformation is `log10(x + 1)`, matching the other lipidomics workflows in this repository.

## What the script does

The script:

- reads the user-defined Excel sheet
- cleans column names with `janitor::clean_names()`
- checks for `infection` and `condition`
- identifies individual lipid columns
- converts lipid columns to numeric
- applies `log10(x + 1)` transformation when `log_transform <- TRUE`
- removes lipid columns with only missing values
- removes lipid columns with zero variance
- keeps complete cases for the lipid matrix and metadata
- runs PERMANOVA using `infection * condition`
- saves the PERMANOVA result
- tests dispersion across combined `infection:condition` groups
- saves the dispersion test result

## PERMANOVA model

The PERMANOVA model is:

```r
resp_mat_cc ~ infection * condition
```

using:

```r
vegan::adonis2(
  method = "euclidean",
  permutations = 999,
  by = "term"
)
```

This tests:

- the infection main effect
- the condition main effect
- the infection x condition interaction

## Dispersion check

The script also runs:

```r
vegan::betadisper()
vegan::permutest()
```

using combined `infection:condition` groups.

This checks whether groups differ in within-group dispersion. This is important because PERMANOVA can be affected by dispersion differences as well as location differences.

## Output files

Outputs are written to `output_dir`.

The script writes:

```text
sms12_individual_lipids_analysis_summary.txt
sms12_individual_lipids_permanova.txt
sms12_individual_lipids_dispersion_test.txt
```

## Main R objects created

The script creates:

- `df_sms12_ind`: cleaned and transformed individual lipid data
- `lipid_cols`: individual lipid measurement columns kept after QC
- `resp_mat_cc`: complete-case lipid response matrix
- `df_cc`: complete-case metadata and transformed lipid data
- `sms12_ind_adonis`: PERMANOVA result
- `dist_ind`: Euclidean distance object
- `disp_ind`: betadisper model
- `disp_perm_ind`: permutation test for dispersion

## Notes

- This script is designed for individual lipid species, not total lipid classes.
- Classical MANOVA is avoided here because high-dimensional individual lipid matrices can exceed the available sample-size and rank constraints.
- PERMANOVA is run on the log-transformed lipid matrix.
- Dispersion is tested on the same complete-case lipid matrix used for PERMANOVA.
- Interpret PERMANOVA together with the dispersion test.
