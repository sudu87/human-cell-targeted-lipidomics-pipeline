# PCA and PERMANOVA Hela vs HUVEC analysis

This document explains how to use the script `scripts/pca_permanova_hela_huvec_analysis.R`.

## Instructions

- Place your input Excel file in your project folder.
- Open `scripts/pca_permanova_hela_huvec_analysis.R`.
- Replace the placeholder input file path with the path to your Excel file:

```r
input_file <- "path/to/your/input_file.xlsx"
```

- Replace the placeholder output directory with the folder where you want results to be saved:

```r
out_dir <- "path/to/your/output_directory/"
```

- Make sure your input file is an `.xlsx` file.
- Make sure the workbook contains metadata columns named `infection` and `cell_line`.
- Keep rows as individual samples and columns as metadata plus lipid measurements.
- Make sure lipid measurement columns contain numeric values.
- Make sure infection and cell-line labels are written consistently.
- Run the script after loading the required R packages.
- Review the PERMANOVA result, PCA plot, scree plot, and output files written to the output directory.

## Overview

This script performs data preparation, PERMANOVA, PCA, and PCA plotting for a lipidomics dataset containing two metadata variables: infection status and cell line.

The script then:

- reads one Excel file
- cleans column names
- encodes `infection` and `cell_line` as categorical factors
- identifies lipid measurement columns
- applies `log10(x + 1)` transformation
- removes samples with incomplete metadata or lipid values
- removes lipids with zero variance
- runs PERMANOVA for `infection`, `cell_line`, and their interaction
- runs PCA
- saves PCA variance explained
- creates a PCA PC1/PC2 plot
- creates a PCA scree plot

## Input requirements

Your Excel file must:

- be an `.xlsx` file
- contain a metadata column named `infection`
- contain a metadata column named `cell_line`
- contain lipid measurements in the remaining columns
- have rows corresponding to individual samples
- have columns corresponding to metadata and lipid measurements

## Specific lines you may need to edit

Update the input file:

```r
input_file <- "path/to/your/input_file.xlsx"
```

Update the output directory:

```r
out_dir <- "path/to/your/output_directory/"
```

Adjust the distance method if needed:

```r
METHOD <- "euclidean"
```

If needed, update the PCA point shapes here:

```r
scale_shape_manual(values = c(21, 24, 22, 23, 25))
```

## What the script does

### Data preparation

- reads the Excel file
- cleans column names
- converts `infection` and `cell_line` to factors
- identifies lipid columns
- applies `log10(x + 1)` transformation to lipid values

### Filtering

- removes samples with missing lipid values or missing metadata
- removes lipid columns with zero variance
- creates:

  - `X` for the filtered lipid matrix
  - `meta` for the filtered sample metadata

### PERMANOVA

- runs PERMANOVA with `adonis2`
- tests the terms `infection`, `cell_line`, and `infection:cell_line`
- uses 999 permutations
- writes the PERMANOVA result to a text file

### PCA

- runs PCA on the filtered lipid matrix
- centers and scales lipid values for PCA
- calculates variance explained by each principal component
- writes variance explained to a `.csv` file
- creates a PC1/PC2 PCA plot
- creates a scree plot

## Output

The script writes the following output files to the output directory:

- `permanova_results.txt`
- `pca_variance_explained.csv`
- `pca_pc1_pc2.png`
- `pca_scree_plot.png`

It also creates the following R objects in the session:

- `df`
- `lipid_cols`
- `X`
- `meta`
- `permanova_fit`
- `pca_fit`
- `pca_var_df`
- `scores`
- `p_pca_2`
- `p_scree`

## Notes

- The script assumes `infection` and `cell_line` are metadata columns.
- All other columns are treated as lipid measurements.
- Non-numeric lipid values may be converted to `NA`.
- Samples with missing lipid values are removed before PCA and PERMANOVA.
- PCA is centered and scaled using `prcomp(..., center = TRUE, scale. = TRUE)`.
- PERMANOVA is run on the log-transformed lipid matrix.
