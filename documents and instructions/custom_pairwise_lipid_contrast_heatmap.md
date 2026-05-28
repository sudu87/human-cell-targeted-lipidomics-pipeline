# Custom pairwise lipid contrast heatmap

This document explains how to use `scripts/custom_pairwise_lipid_contrast_heatmap.R`.

The script creates a custom ordered, FDR-filtered heatmap for selected pairwise lipid contrasts. It can be used for either total lipid classes or individual lipid species, as long as the input file contains one row per lipid and contrast.

## Required R packages

The script uses:

```r
readxl
dplyr
ggplot2
```

Install missing packages before running the script.

## User-defined paths

Open `scripts/custom_pairwise_lipid_contrast_heatmap.R` and set:

```r
input_file <- "path/to/your/pairwise_results.xlsx"
sheet_name <- NULL
output_dir <- "path/to/your/output_directory"
```

Use `sheet_name <- NULL` to read the first worksheet, or set `sheet_name` to a specific Excel sheet name.

## Input data requirements

The input file should contain pairwise lipid contrast results with columns for:

- lipid name
- contrast name
- contrast estimate
- FDR-adjusted p-value
- infection status

By default, the script expects these column names:

```r
lipid_col <- "lipid"
contrast_col <- "contrast"
estimate_col <- "estimate"
p_value_col <- "p_fdr"
infection_col <- "infection"
```

If your file uses different names, update these variables before running the script.

## Choosing pairwise comparisons

The pairwise comparisons shown in the heatmap are chosen by:

```r
keep_contrasts <- c(
  "d_sms12_DKO - h_sms12_WT",
  "a_sms1_AHT - e_sms1_no_AHT",
  "e_sms1_no_AHT - h_sms12_WT",
  "a_sms1_AHT - h_sms12_WT",
  "b_sms2_AHT - f_sms2_no_AHT",
  "f_sms2_no_AHT - h_sms12_WT",
  "b_sms2_AHT - h_sms12_WT",
  "(c_sms2-M64R_AHT) - (g_sms2-M64R_no_AHT)",
  "(g_sms2-M64R_no_AHT) - h_sms12_WT",
  "(c_sms2-M64R_AHT) - h_sms12_WT"
)
```

This vector does two things:

- chooses which pairwise comparisons are kept
- sets the left-to-right x-axis order of those comparisons

Only contrasts listed in `keep_contrasts` are plotted.

## Choosing lipid order

The lipid rows shown in the heatmap are chosen by:

```r
custom_lipid_order <- c(
  "dh_sph", "sph", "s1p",
  "dh_cer16_0", "dh_cer18_0",
  "cer16_0", "cer18_0",
  "dh_sm16_0", "dh_sm18_0",
  "sm16_0", "sm18_0",
  "hex_cer16_0",
  "lac_cer16_0"
)
```

This vector does two things:

- chooses which lipids are kept when `keep_only_custom_lipids <- TRUE`
- sets the top-to-bottom y-axis order

For total lipid classes, use a shorter order such as:

```r
custom_lipid_order <- c(
  "dh_sph",
  "sph",
  "s1p",
  "dh_cer_total",
  "cer_total",
  "dh_sm_total",
  "sm_total",
  "hex_cer_total",
  "lac_cer_total"
)
```

For individual lipid species, use a longer species-level order such as:

```r
custom_lipid_order <- c(
  "dh_sph", "sph", "s1p",
  "dh_cer16_0", "dh_cer18_0", "dh_cer20_0",
  "cer16_0", "cer18_0", "cer20_0",
  "dh_sm16_0", "dh_sm18_0", "dh_sm20_0",
  "sm16_0", "sm18_0", "sm20_0"
)
```

Set:

```r
keep_only_custom_lipids <- FALSE
```

to append any additional lipids after the custom-ordered lipids.

## Estimate scale

The script converts the contrast estimate to log2 fold-change using:

```r
estimate_scale <- "log10"
```

This is appropriate when the pairwise model was run on `log10(x + 1)` lipid values. Other supported options are:

```r
estimate_scale <- "log2"
estimate_scale <- "natural_log"
estimate_scale <- "fold_change"
```

## What the script does

The script:

- reads the pairwise contrast Excel file
- standardizes the configured input columns
- keeps only `keep_contrasts`
- keeps only lipids in `custom_lipid_order` when `keep_only_custom_lipids <- TRUE`
- converts estimates to log2 fold-change
- keeps significant values where the configured p-value column is below `alpha`
- sets non-significant values to `NA`
- plots significant log2 fold-change values in a heatmap
- shows non-significant or missing cells as grey
- facets the heatmap by infection status
- writes the plot data and missing contrast/lipid diagnostics
- saves the heatmap as PDF and PNG

## Output files

Outputs are written to `output_dir`.

The script writes:

```text
custom_pairwise_lipid_contrast_heatmap_plot_data.csv
custom_pairwise_lipid_contrast_heatmap_missing_contrasts.txt
custom_pairwise_lipid_contrast_heatmap_missing_lipids.txt
custom_pairwise_lipid_contrast_heatmap.pdf
custom_pairwise_lipid_contrast_heatmap.png
```

## Main R objects created

The script creates:

- `df_heat`: raw pairwise contrast table
- `df_hm`: filtered heatmap input table
- `plot_df`: transformed plot data with log2 fold-change and significance filter
- `p_heatmap`: ggplot heatmap object

## Notes

- This script can be used for total lipid classes or individual lipid species.
- The key customization is `keep_contrasts`, which chooses and orders the pairwise comparisons.
- The second key customization is `custom_lipid_order`, which chooses and orders lipid rows.
- Non-significant cells are grey because `log2FC_sig` is set to `NA` when `p_value >= alpha`.
- X-axis labels are hidden by default. Set `show_x_axis_labels <- TRUE` to display contrast names.
