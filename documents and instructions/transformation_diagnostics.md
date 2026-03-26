
# Transformation diagnostics

This document explains how to use the script `scripts/transformation_diagnostics.R`.

## Instructions

- Make sure the main analysis script has already created the objects `df`, `lipid_cols`, and `out_dir`.
- Open `scripts/transformation_diagnostics.R`.
- Replace the placeholder output directory in the following line if needed:

```r
out_dir <- "path/to/your/output_directory/"
```

* Make sure `df` contains the lipid measurement data.
* Make sure `lipid_cols` contains the names of the lipid measurement columns.
* Make sure your dataset contains a column named `infection`, because it is used for the ANOVA residual check.
* Make sure lipid measurement columns contain numeric values.
* Run the script after loading the required R packages.
* Check the generated diagnostic `.png` files in the diagnostics output folder.

## Overview

This script generates diagnostic plots to inspect the distribution of lipid values before and after log10 transformation, and to assess pooled ANOVA residual normality.

## Input requirements

The script expects the following objects to be available in the R session:

* `df`: data frame containing sample metadata and lipid measurements
* `lipid_cols`: character vector containing the lipid measurement column names
* `out_dir`: output directory for analysis results

Your data must also contain:

* a column named `infection`
* numeric lipid measurement columns

## Specific lines you may need to edit

Update the output directory if needed:

```r
out_dir <- "path/to/your/output_directory/"
```

If this script is run independently, make sure `df` and `lipid_cols` are defined before running it.

## What the script does

The script:

* creates a `diagnostics` subfolder inside the main output directory
* extracts all lipid values into one vector
* plots a histogram of raw lipid values
* plots a histogram of log10-transformed lipid values
* plots density curves for raw and log10-transformed values
* calculates pooled residuals from one-way ANOVA models using `infection` as predictor
* generates a QQ-plot of pooled ANOVA residuals
* saves all plots as `.png` files

## Output

The script writes the following files to the diagnostics output directory:

- `hist_raw_values.png`
- `hist_log10_values.png`
- `density_raw_vs_log.png`
- `qqplot_residuals.png`

## Example plots

### Raw value distribution
![Raw value distribution](images/hist_raw_values.png)

### Log10-transformed value distribution
![Log10-transformed value distribution](images/hist_log10_values.png)

### Density plots
![Density plots](images/density_raw_vs_log.png)

### QQ-plot of pooled residuals
![QQ-plot of pooled residuals](images/qqplot_residuals.png)

## Notes

* The script assumes `infection` is available in `df`.
* The QQ-plot uses pooled residuals across all lipid-wise ANOVA models.
* The script does not return a main summary object.
* Make sure the required packages and data objects are available before running the script.

