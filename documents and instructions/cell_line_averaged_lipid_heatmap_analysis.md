# Averaged lipid heatmap analysis

This document explains how to use the script `scripts/averaged_lipid_heatmap_analysis.R`.

## Instructions

- Place your input Excel file in your project folder.
- Open `scripts/averaged_lipid_heatmap_analysis.R`.
- Replace the placeholder input file path with the path to your Excel file:

```r
input_file <- "path/to/your/input_file.xlsx"
```

- Replace the placeholder sheet name with the worksheet you want to analyze:

```r
sheet_name <- "your_sheet_name"
```

- Replace the placeholder output directory with the folder where you want the results to be saved:

```r
out_dir <- "path/to/your/output_directory/"
```

- Make sure the Excel file contains metadata columns named `infection` and `cell_line`.
- Keep rows as individual samples and columns as metadata plus lipid measurements.
- Update the factor order if your labels are different:

```r
cell_line_order <- c("cell_line_1", "cell_line_2")
infection_order <- c("uninfected", "infected")
```

- Make sure lipid measurement columns contain numeric values.
- Run the script after loading the required R packages.
- Review the averaged lipid table and heatmap image saved in the output directory.

## Overview

This script creates a heatmap from lipidomics data after averaging samples by cell line and infection group.

The script then:

- reads one Excel worksheet
- cleans column names
- defines metadata and lipid columns
- converts `infection` and `cell_line` to factors
- converts lipid columns to numeric
- applies `log10(x + 1)` transformation
- averages lipid values by `cell_line` and `infection`
- z-scores each lipid across averaged groups
- creates a heatmap using `pheatmap`
- saves the averaged lipid table and heatmap image

## Input requirements

Your Excel file must:

- be provided as an `.xlsx` file
- contain a worksheet with the lipidomics data
- contain metadata columns named `infection` and `cell_line`
- contain lipid measurements in the remaining columns
- have rows corresponding to individual samples
- have columns corresponding to metadata and lipid measurements

## Specific lines you may need to edit

Update the input file:

```r
input_file <- "path/to/your/input_file.xlsx"
```

Update the worksheet name:

```r
sheet_name <- "your_sheet_name"
```

Update the output directory:

```r
out_dir <- "path/to/your/output_directory/"
```

Update the factor order:

```r
cell_line_order <- c("cell_line_1", "cell_line_2")
infection_order <- c("uninfected", "infected")
```

If needed, update the heatmap colours:

```r
pal <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(101)
```

If needed, add fixed annotation colours:

```r
annotation_colours <- list(
  infection = c("uninfected" = "#009E73", "infected" = "#CC79A7"),
  cell_line = c("cell_line_1" = "#E69F00", "cell_line_2" = "#999999")
)
```

## What the script does

### Data preparation

- reads the Excel worksheet
- cleans column names
- identifies metadata columns
- treats all remaining columns as lipid measurements
- converts lipid measurements to numeric
- applies `log10(x + 1)` transformation

### Averaging

- groups samples by `cell_line` and `infection`
- calculates the mean value for each lipid within each group
- creates a sample ID from the combined cell-line and infection labels
- saves the averaged values as a `.csv` file

### Heatmap

- builds a lipid matrix from the averaged values
- removes lipids with zero variance
- z-scores each lipid across averaged groups
- transposes the matrix so samples appear as heatmap columns
- adds sample annotation for infection and cell line
- creates and saves the heatmap

## Output

The script writes the following files to the output directory:

- `averaged_lipid_values.csv`
- `averaged_lipid_heatmap.png`

It also creates the following R objects in the session:

- `df`
- `lipid_cols`
- `df_avg`
- `mat`
- `mat_t`
- `annotation_col`
- `heatmap_obj`

## Notes

- The script averages lipid values after log transformation.
- Z-scoring is performed after averaging.
- Lipids with zero variance are removed before z-scoring.
- The names in `annotation_colours` must match the factor labels used in your data.
- If you do not need fixed annotation colours, keep `annotation_colours <- NULL`.
