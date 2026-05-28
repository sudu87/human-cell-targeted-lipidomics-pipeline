# SMS1/2 selected lipid summary bar plots

This document explains how to use `scripts/sms12_lipid_summary_barplots.R`.

The script creates grouped bar plots for selected SMS1/2 lipid totals. It reads a raw Excel file, summarises each selected lipid by `condition` and `infection`, and writes both individual lipid plots and one combined faceted plot.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
tidyr
ggplot2
```

Install missing packages before running the script.

## User-defined paths

Open `scripts/sms12_lipid_summary_barplots.R` and set:

```r
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"
```

The input file should be the raw Excel workbook containing the lipid measurements.

## Input data requirements

After `janitor::clean_names()` is applied, the Excel sheet must contain:

- one row per sample
- a column named `infection`
- a column named `condition`
- numeric lipid measurement columns

By default, the script plots:

```r
lipids_interest <- c(
  "cer_total",
  "sm_total",
  "hex_cer_total",
  "lac_cer_total",
  "dh_sm_total",
  "sph"
)
```

Edit `lipids_interest` if your input file uses different lipid columns. If a selected lipid is missing, the script warns and skips it. If none of the selected lipids are found, the script stops.

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
- reshapes selected lipids from wide format to long format
- calculates mean, standard deviation, sample size, and standard error for each `condition x infection x lipid` group
- writes the summary table as a CSV file
- creates one bar plot per selected lipid
- creates one combined faceted plot with a shared legend

The plotted values are raw lipid values, not log-transformed values.

## Output files

Outputs are written to `output_dir`.

The script writes:

```text
sms12_lipid_summary_statistics.csv
sms12_selected_lipid_barplots.pdf
sms12_selected_lipid_barplots.png
<lipid_name>_barplot.pdf
<lipid_name>_barplot.png
```

For example:

```text
cer_total_barplot.pdf
sm_total_barplot.png
```

## Main R objects created

The script creates:

- `df_sms12_raw`: cleaned raw input data
- `lipids_present`: selected lipid columns found in the input file
- `df_long`: long-format lipid data
- `sum_df`: summary statistics by condition, infection, and lipid
- `plots_by_lipid`: named list of individual ggplot objects
- `p_combined`: combined faceted ggplot object

## Notes

- The original plotting idea used separate plot objects and `gridExtra::grid.arrange()`.
- This standalone script uses `facet_wrap()` for the combined plot, which avoids requiring extra packages such as `cowplot` or `gridExtra`.
- The y-axis starts at zero for the bar plots.
- Error bars show mean ± standard error.
- If your raw file has lipid columns with different cleaned names, update `lipids_interest` and `lipid_labels`.
