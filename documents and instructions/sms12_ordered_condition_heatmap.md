# SMS1/2 ordered-condition lipid heatmap

This document explains how to use `scripts/sms12_ordered_condition_heatmap.R`.

The script creates SMS1/2 lipid heatmaps for user-selected infection subsets. It is a special-case heatmap workflow where the user defines the exact condition order and optional display labels for the heatmap columns.

## Required R packages

The script uses:

```r
readxl
janitor
dplyr
pheatmap
```

Install missing packages before running the script.

## User-defined paths

Open `scripts/sms12_ordered_condition_heatmap.R` and set:

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
- lipid measurement columns in the remaining columns

The script defines lipid columns as all columns except:

```r
meta_cols <- c("infection", "condition")
```

## Infection subsets

The infection groups to plot are controlled by:

```r
infection_groups_to_plot <- c("yes", "no")
```

Use:

```r
infection_groups_to_plot <- c("yes")
```

to create only the infected heatmap, or:

```r
infection_groups_to_plot <- c("no")
```

to create only the uninfected heatmap.

The original exploratory code created `df_yes` and `df_no`, but then transformed `df_no` while plotting `df_yes`. This standalone script fixes that bug by filtering one infection group at a time and applying numeric conversion and `log10(x + 1)` transformation to the same subset that is plotted.

## Condition order and labels

The heatmap column order is controlled by:

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

Conditions found in the data but not listed in `condition_order` are appended after the ordered conditions.

Column display labels are controlled by:

```r
condition_labels <- c(
  "sms12_WT" = "HeLa WT (-AHT)",
  "sms12_DKO" = "SMS1/2 DKO (-AHT)",
  "sms1_16_no_AHT" = "SMS1 comp (-AHT)",
  "sms1_16_AHT" = "SMS1 comp (+AHT)",
  "sms2_no_AHT" = "SMS2 comp (-AHT)",
  "sms2_AHT" = "SMS2 comp (+AHT)",
  "sms2-M64R_no_AHT" = "SMS2-M64R (-AHT)",
  "sms2-M64R_AHT" = "SMS2-M64R (+AHT)"
)
```

Any condition without a label keeps its original condition name.

## What the script does

For each infection group in `infection_groups_to_plot`, the script:

- filters the raw data to that infection group
- converts lipid columns to numeric
- applies `log10(x + 1)` transformation when `log_transform <- TRUE`
- averages lipid values across replicate rows within each condition
- writes the averaged values to a CSV file
- builds a matrix with lipids as rows and conditions as columns
- applies the user-defined condition order
- applies the user-defined condition labels
- removes lipids with no finite values
- removes zero-variance lipids before row scaling
- draws a heatmap with row-wise z-scoring
- saves PNG and PDF heatmap files

## Output files

Outputs are written to `output_dir`.

For each infection value, the script writes:

```text
sms12_ordered_condition_heatmap_<infection>_averaged_values.csv
sms12_ordered_condition_heatmap_<infection>.png
sms12_ordered_condition_heatmap_<infection>.pdf
```

For example, with `infection_groups_to_plot <- c("yes", "no")`, the script writes:

```text
sms12_ordered_condition_heatmap_yes_averaged_values.csv
sms12_ordered_condition_heatmap_yes.png
sms12_ordered_condition_heatmap_yes.pdf
sms12_ordered_condition_heatmap_no_averaged_values.csv
sms12_ordered_condition_heatmap_no.png
sms12_ordered_condition_heatmap_no.pdf
```

## Main R objects created

The script creates:

- `df`: cleaned raw input data
- `lipid_cols`: lipid measurement columns
- `make_heatmap_for_infection`: helper function that creates one infection-specific heatmap
- `heatmap_results`: named list containing outputs for each infection group

Each entry in `heatmap_results` contains:

- `infection`: infection value used for that heatmap
- `averaged_values`: averaged condition-level lipid data
- `matrix`: heatmap matrix after ordering and labelling
- `heatmap_png`: pheatmap object for the PNG output
- `heatmap_pdf`: pheatmap object for the PDF output

## Notes

- The plotted matrix is row-scaled by `pheatmap(scale = "row")`, so each lipid is z-scored across conditions.
- Rows and columns are not clustered by default.
- The heatmap colour scale is fixed from -2 to 2.
- The script corrects the original subset bug by transforming the same infection subset that is plotted.
- The `sms2_AHT` display label is set to `SMS2 comp (+AHT)`.
