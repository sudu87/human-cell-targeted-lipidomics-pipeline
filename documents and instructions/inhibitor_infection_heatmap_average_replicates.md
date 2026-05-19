

# Heatmap from averaged replicates

This document explains how to use the script `scripts/inhibitor_infection_heatmap_average_replicates.R`.

## Instructions

- Place your input Excel file in your project folder or note its full file path.
- Open `scripts/heatmap_average_replicates.R`.
- Replace the placeholder path in the following line with the path to your own Excel file:

```r
file_path <- "path/to/your/input_file.xlsx"
```

* Make sure your Excel file is an `.xlsx` file.
* Make sure your Excel file contains the metadata columns `infection` and `inhibitor`.
* Arrange all lipid species measurements in the remaining columns.
* Keep replicate samples as separate rows, because the script averages them automatically.
* Check that the `infection` column contains the exact labels `Uninfected` and `infected`.
* Check that the `inhibitor` column contains the expected labels `untreated`, `DMSO`, `AKS466`, `HPA-12`, `Desipramin`, `Myriocin`, and `ARC39`.
* Edit the factor levels in the script only if your dataset uses different condition names or a different sample order.
* Make sure the lipid measurement columns contain numeric values, because non-numeric entries may be converted to `NA`.
* Run the script after loading the required R packages.
* Use the resulting object `heatmap_object` for later plotting, exporting, or combining with other figures.

## Overview

This script generates a heatmap from lipid measurements after averaging replicates by `inhibitor` and `infection`. The heatmap is stored as an R object for later use.

## Input requirements

Your Excel file must:

* be provided as an `.xlsx` file
* contain the metadata columns `infection` and `inhibitor`
* contain lipid species measurements in all remaining columns
* keep replicate samples as separate rows
* use the expected labels in the metadata columns

Expected values:

* `infection`: `Uninfected`, `infected`
* `inhibitor`: `untreated`, `DMSO`, `AKS466`, `HPA-12`, `Desipramin`, `Myriocin`, `ARC39`

## Specific lines you may need to edit

Update this line with the path to your own Excel file:

```r
file_path <- "path/to/your/input_file.xlsx"
```

Edit these lines only if your dataset uses different condition names or a different order:

```r
inhibitor = factor(inhibitor, levels = c("untreated","DMSO","AKS466","HPA-12","Desipramin","Myriocin","ARC39")),
infection = factor(infection, levels = c("Uninfected","infected"))
```

## What the script does

The script:

* reads the Excel file
* standardizes column names
* identifies metadata and lipid columns
* converts lipid columns to numeric values
* applies `log10(x + 1)` transformation before averaging
* averages replicate rows by `inhibitor` and `infection`
* enforces a consistent sample order
* z-scores lipid values across samples
* creates a heatmap with fixed colors and breaks
* stores the heatmap in the object `heatmap_object`

## Output

The final heatmap is stored in:

```r
heatmap_object
```

This allows you to reuse the heatmap later for plotting, exporting, or combining with other figures.

## Notes

* Non-numeric values in lipid columns may be converted to `NA`.
* Zero values are allowed because the script uses `log10(x + 1)`.
* The script does not save the heatmap to a file automatically.
* Make sure the required packages are loaded before running the script.

````


