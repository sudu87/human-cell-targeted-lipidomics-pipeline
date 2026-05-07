# MANOVA followed by Per-Lipid ANOVA and Post-hoc analysis

This document explains how to use the script `scripts/manova_per_lipid_anova_analysis.R`.

## Background 
## When MANOVA is appropriate

MANOVA is appropriate when:

- Multiple lipid species are measured from the same experimental unit.
- The lipid variables are biologically correlated.
- The goal is to test an overall lipidome-level shift.
- The samples are independent biological replicates.
- The number of response variables is reasonable relative to the residual degrees of freedom.

This is often conceptually appropriate for sphingolipid or lipidomics data because lipid species within the same pathway are typically correlated and are not independent biological endpoints.


## When MANOVA is not recommended

Avoid or interpret MANOVA cautiously when:

- The main goal is only to interpret each lipid individually.
- The number of lipids is larger than, or close to, the residual degrees of freedom.
- The design is very unbalanced or has missing infection-by-inhibitor cells.
- There are very few biological replicates per group.
- Many lipid variables are nearly perfectly correlated, causing rank-deficient residual covariance matrices.
- Distributions are highly skewed even after transformation.

For high-dimensional lipidomics data, consider reducing the response matrix before classical MANOVA, for example by filtering lipid species, using pathway-level summaries, or testing principal components. A permutation-based multivariate method such as PERMANOVA may also be appropriate for exploratory global lipidome testing, but that is a different null hypothesis and should be reported separately.


## Why Pillai's trace?

Pillai's trace is one of the standard MANOVA test statistics, alongside Wilks' lambda, Hotelling–Lawley trace, and Roy's largest root. The R `summary.manova()` method supports all four and uses Pillai as the first/default option. Pillai's trace is often preferred as a robust default when biological data may show moderate deviations from ideal MANOVA assumptions.

Interpretation:

- Higher Pillai's trace indicates stronger multivariate group separation.
- A small p-value, for example `p < 0.05`, indicates evidence against the null hypothesis that the multivariate lipid mean vectors are equal for that model term.
- A significant MANOVA result does **not** identify which individual lipid differs; follow-up per-lipid ANOVA and post-hoc contrasts are required.

## Instructions

- Place your input Excel file in your chosen input location.
- Open `scripts/total_lipids_manova_simple.R`.
- Replace the placeholder input file path with your own Excel file:

```r
input_file <- "path/to/input_file.xlsx"
```

- Replace the placeholder output directory with the folder where results should be saved:

```r
out_dir <- "path/to/output_directory"
```

- Make sure the Excel file contains metadata columns named `infection` and `inhibitor`.
- Keep rows as individual samples and columns as metadata plus lipid measurements.
- Make sure lipid measurement columns contain numeric values.
- Run the script after loading or installing the required R packages.

## Overview

This script performs a simple MANOVA workflow for total lipid data.

The script:

- reads one Excel file
- cleans column names
- encodes `infection` and `inhibitor` as factors
- identifies lipid columns
- applies `log10(x + 1)` transformation
- creates a lipid response matrix
- runs MANOVA using `infection * inhibitor`
- reports Pillai's trace
- reports per-lipid ANOVA results using `summary.aov`
- plots residuals with a normal density curve
- runs posthoc infection contrasts within each inhibitor
- runs posthoc inhibitor contrasts within each infection group
- applies BH/FDR correction in the posthoc tests
- writes posthoc results to Excel files

## Input requirements

Your Excel file must:

- be an `.xlsx` file
- contain one row per sample
- contain a column named `infection`
- contain a column named `inhibitor`
- contain lipid measurements in the remaining columns
- use consistent infection labels
- use consistent inhibitor labels

## Specific lines you may need to edit

Update the input file:

```r
input_file <- "path/to/input_file.xlsx"
```

Update the output directory:

```r
out_dir <- "path/to/output_directory"
```

The script saves output files using these generic names:

```r
out_infection <- file.path(out_dir, "posthoc_infection_contrasts_by_inhibitor.xlsx")
out_inhibitor <- file.path(out_dir, "posthoc_inhibitor_contrasts_by_infection.xlsx")
```

## What the script does

### MANOVA

The script creates a lipid matrix and fits:

```r
manova_fit <- manova(resp_mat ~ infection * inhibitor, data = df)
```

The MANOVA result is printed using Pillai's trace:

```r
summary(manova_fit, test = "Pillai")
```

Per-lipid ANOVA results are printed using:

```r
summary.aov(manova_fit)
```

### Posthoc analysis 1

The first posthoc analysis compares infection groups within each inhibitor:

```r
em <- emmeans(fit, ~ infection | inhibitor)
pw <- pairs(em, adjust = "BH")
```

### Posthoc analysis 2

The second posthoc analysis compares inhibitors within each infection group:

```r
em <- emmeans(fit, ~ inhibitor | infection)
pw <- pairs(em, adjust = "BH")
```

## Output

The script writes two Excel files:

- `posthoc_infection_contrasts_by_inhibitor.xlsx`
- `posthoc_inhibitor_contrasts_by_infection.xlsx`

It also prints the MANOVA and per-lipid ANOVA results in the R console.


## Notes

- The posthoc tests use `adjust = "BH"`, so the p-values are BH/FDR corrected.
- Pillai's trace is used for the main MANOVA test.
- The script uses `log10(x + 1)` transformation before MANOVA and ANOVA.
- `plotNormalDensity()` is used only as a visual residual check.
- The script uses the `infection` and `inhibitor` labels already present in the Excel file.
