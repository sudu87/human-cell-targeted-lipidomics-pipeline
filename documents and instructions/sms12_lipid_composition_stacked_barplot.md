# SMS1/2 lipid composition stacked bar plot

This document explains how to use `scripts/sms12_lipid_composition_stacked_barplot.R`.

The script creates 100% stacked bar plots showing the relative composition of selected sphingolipid classes across SMS1/2 conditions. It uses the same type of raw input data as the SMS1/2 lipid summary barplot workflow, but plots proportions rather than absolute abundance.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
tidyr
ggplot2
scales
```

Install missing packages before running the script.

## User-defined paths

Open `scripts/sms12_lipid_composition_stacked_barplot.R` and set:

```r
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"
```

The input file should be the raw Excel workbook containing total lipid measurements.

## Input data requirements

After `janitor::clean_names()` is applied, the Excel sheet must contain:

- one row per sample
- a column named `infection`
- a column named `condition`
- lipid measurement columns

The default lipid classes are:

```r
lipids_interest <- c(
  "cer_total",
  "sm_total",
  "hex_cer_total",
  "lac_cer_total",
  "dh_sm_total"
)
```

Edit `lipids_interest`, `lipid_labels`, and `lipid_palette` if your input file uses different lipid names or if you want different labels/colours.

## Condition and infection order

The x-axis order is controlled by:

```r
condition_order <- c(
  "sms12_WT",
  "sms12_DKO",
  "sms1_16_no_AHT",
  "sms1_16_AHT",
  "sms2_no_AHT",
  "sms2_AHT",
  "sms2-M64R_no_AHT",
  "sms2-M64R_AHT"
)
```

Set `condition_order <- NULL` to use the condition order found in the input file.

The infection groups are controlled by:

```r
infection_order <- c("no", "yes")
infection_labels <- c("no" = "Uninfected", "yes" = "Infected")
```

Update these if your input file uses different infection labels.

## What the script does

The script:

- reads the user-defined Excel sheet
- cleans column names with `janitor::clean_names()`
- checks for `infection` and `condition`
- converts selected lipid columns to numeric
- averages selected lipid values by `condition x infection`
- reshapes the averaged data from wide format to long format
- calculates each selected lipid as a proportion of the selected lipid total for that `condition x infection` group
- creates a 100% stacked bar plot with one facet per infection group
- writes the averaged and long-format data tables
- saves the stacked bar plot as PDF and PNG

The plotted values are proportions calculated from raw lipid abundances. The script does not log-transform the values because proportions should reflect the original lipid composition.

## Output files

Outputs are written to `output_dir`.

The script writes:

```text
sms12_lipid_composition_averaged_values.csv
sms12_lipid_composition_long_values.csv
sms12_lipid_composition_stacked_barplot.pdf
sms12_lipid_composition_stacked_barplot.png
```

## Main R objects created

The script creates:

- `df_sms12_raw`: cleaned raw input data
- `lipids_present`: selected lipid columns found in the input file
- `df_mean_inf`: mean lipid abundance by condition and infection
- `df_long_mean_inf`: long-format data with lipid proportions
- `p_stack_inf`: ggplot stacked bar plot object

## Notes

- This plot answers a composition question: what proportion of the selected sphingolipid pool is represented by each lipid class?
- It is different from the grouped barplot workflow, which plots absolute mean abundance for each lipid.
- Missing selected lipids are skipped with a warning.
- The y-axis uses percentages through `scales::percent_format()`.
