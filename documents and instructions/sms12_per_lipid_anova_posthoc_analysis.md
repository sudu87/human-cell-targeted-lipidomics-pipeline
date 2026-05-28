# Per-lipid ANOVA and posthoc analysis

This document explains how to use `scripts/sms12_per_lipid_anova_posthoc_analysis.R`.

The script performs per-lipid linear models, ANOVA tables, posthoc contrasts, and FDR correction for lipidomics data.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
purrr
tidyr
tibble
emmeans
writexl
```

Install missing packages before running the script.

## Configure input and output paths

Before running the script, edit the configuration block near the top of the script:

```r
input_file <- "path/to/input_lipidomics_data.xlsx"
sheet_name <- "sheet_name"
out_dir <- "analysis_outputs/statistics"
```

Replace `input_file` with the path to your Excel workbook, `sheet_name` with the sheet to analyze, and `out_dir` with the folder where output files should be saved.

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
- fits one linear model per lipid using `infection * condition`
- exports per-lipid ANOVA results
- computes pairwise condition contrasts within each infection group
- computes pairwise infection contrasts within each condition
- standardizes the posthoc p-value column as `p.value`
- applies Benjamini-Hochberg FDR correction across lipids for the same contrast within each stratum
- exports posthoc results to Excel

## Output files

Outputs are written to the folder defined by:

```r
out_dir
```

The script writes:

```text
per_lipid_anova.xlsx
per_lipid_posthoc.xlsx
```

## Main R objects created

The script creates:

- `sms_test`: cleaned and log-transformed input data
- `lipid_cols`: lipid measurement columns
- `anova_sms_test`: combined per-lipid ANOVA table
- `res`: unadjusted posthoc contrast results
- `cond_tbl`: condition contrasts within infection groups with FDR correction
- `inf_tbl`: infection contrasts within condition groups with FDR correction
- `posthoc_sms_test`: combined FDR-corrected posthoc results

## Notes

- The per-lipid model formula is `lipid ~ infection * condition`.
- Posthoc contrasts use `adjust = "none"` in `emmeans`; FDR correction is then applied manually across lipids for the same contrast within each stratum.
- Condition contrasts are grouped by `infection` and `contrast` before FDR correction.
- Infection contrasts are grouped by `condition` and `contrast` before FDR correction.
