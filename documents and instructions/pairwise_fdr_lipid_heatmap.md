
# Pairwise comparison and FDR heatmaps for lipids

This document explains how to use the script `scripts/pairwise_fdr_lipid_heatmap.R`.

## Instructions

- Place your input Excel file in your project folder or note its full file path.
- Open `scripts/differential_lipid_heatmap.R`.
- Replace the placeholder path in the following line with the path to your own Excel file:

```r
file_path <- "path/to/your/input_file.xlsx"
```

* Make sure your input file is an `.xlsx` file.
* Make sure your dataset contains the columns `lipid`, `estimate`, `p.value`, `infection`, and `contrast`.
* Check the available contrast names using:

```r
dput(unique(df_heat$contrast))
```

* Replace the placeholder lipid names in `lipid_order` with your own lipid order:

```r
lipid_order <- c(
  "lipid1",
  "lipid2"
)
```

* Edit the legend title in the following line so it matches your comparison:

```r
legend_title <- "Treatment\nvs\nUntreated"
```

* Check that the `infection` column contains the exact labels `Uninfected` and `infected`.
* Make sure the `lipid` column contains names that match the entries in your `lipid_order` object.
* Make sure `estimate` and `p.value` contain numeric values.
* Run the script after loading the required R packages.
* Use the resulting object `heatmap_plot` for later plotting, exporting, or combining with other figures.

## Overview

This script generates a heatmap from differential lipid analysis results. It converts `estimate` values to log2 values, masks non-significant values using `p.value < 0.05`, applies a custom lipid order, and stores the final plot as an R object.

## Input requirements

Your Excel file must:

* be provided as an `.xlsx` file
* contain the columns `lipid`, `estimate`, `p.value`, `infection`, and `contrast`
* contain valid numeric values in `estimate` and `p.value`
* use infection labels that match the expected values
* contain lipid names that match the `lipid_order` object

Expected values:

* `infection`: `Uninfected`, `infected`

## Specific lines you may need to edit

Update this line with the path to your own Excel file:

```r
file_path <- "path/to/your/input_file.xlsx"
```

Set your custom lipid order here:

```r
lipid_order <- c(
  "lipid1",
  "lipid2"
)
```

Update the legend title here:

```r
legend_title <- "Treatment\nvs\nUntreated"
```

If needed, edit the infection factor order here:

```r
infection = factor(infection, levels = c("Uninfected", "infected"))
```

## What the script does

The script:

* reads the Excel file
* filters rows with missing `lipid`, `estimate`, `p.value`, `infection`, or `contrast` values
* converts `estimate` values to log2 fold change
* keeps only significant values in the heatmap using `p.value < 0.05`
* sets a custom lipid order
* sorts contrast values
* applies factor ordering to `lipid`, `contrast`, and `infection`
* creates a faceted heatmap by infection status
* stores the heatmap in the object `heatmap_plot`

## Output

The final heatmap is stored in:

```r
heatmap_plot
```

This allows you to reuse the heatmap later for plotting, exporting, or combining with other figures.

## Notes

* Non-significant values are replaced with `NA` in the plotted fill variable.
* These `NA` values are shown in grey in the heatmap.
* The script does not save the heatmap to a file automatically.
* Make sure the required packages are loaded before running the script.

```
```
