# PCA and PERMANOVA inhibitor analysis

This document explains how to use the script `scripts/pca_permanova_inhibitor_analysis.R`.

## Instructions

- Place your input Excel files in one folder.
- Open `scripts/pca_permanova_inhibitor_analysis.R`.
- Replace the placeholder input directory in the following line with the folder containing your Excel files:

```r
data_dir <- "path/to/your/input_directory/"
````

* Replace the placeholder output directory in the following line with the folder where you want results to be saved:

```r
out_dir <- "path/to/your/output_directory/"
```

* Make sure your input files are `.xlsx` files.
* Make sure each Excel file corresponds to one inhibitor condition.
* Make sure each workbook contains a metadata column named `infection`.
* Keep rows as individual samples and columns as metadata plus lipid measurements.
* Update the file name mapping in the script so each input file is assigned to the correct inhibitor label.
* Set `KEEP_ONLY_UNINFECTED <- TRUE` if you want to analyze only uninfected samples.
* Set `KEEP_ONLY_UNINFECTED <- FALSE` if you want to keep all infection groups.
* Make sure the `infection` labels match the expected values used in the script.
* Make sure lipid measurement columns contain numeric values, because non-numeric values may be converted to `NA`.
* Run the script after loading the required R packages.
* Review the PCA plot object and the output files written to the output directory.

## Overview

This script performs data preparation, quality control filtering, PCA, PERMANOVA, and dispersion testing for lipidomics datasets stored across multiple Excel files. Each workbook is treated as one inhibitor condition, the inhibitor label is assigned from the file name, and all files are combined into a single analysis table.

The script then:

* cleans and aligns the imported datasets
* encodes `infection` and `inhibitor` as categorical factors
* optionally subsets to uninfected samples only
* coerces lipid columns to numeric
* applies `log10(x + 1)` transformation
* removes lipids with only missing values
* removes lipids with zero variance
* removes samples with incomplete data
* z-scores lipid values when Euclidean distance is used
* runs PCA
* saves variance explained by each principal component
* creates a PCA plot
* runs PERMANOVA for inhibitor effects
* performs dispersion checks using `betadisper` and `permutest`

## Input requirements

Your Excel files must:

* be provided as `.xlsx` files
* be stored together in one input directory
* contain one workbook per inhibitor condition
* contain a metadata column named `infection`
* contain lipid measurements in the remaining columns
* have rows corresponding to individual samples
* have columns corresponding to metadata and lipid measurements

## Specific lines you may need to edit

Update the input directory:

```r
data_dir <- "path/to/your/input_directory/"
```

Update the output directory:

```r
out_dir <- "path/to/your/output_directory/"
```

Choose whether to keep only uninfected samples:

```r
KEEP_ONLY_UNINFECTED <- TRUE
```

Adjust the analysis method if needed:

```r
METHOD <- "euclidean"
```

Update the file name mapping if your input file names are different:

```r
name_map <- c(
  "AKS466_new" = "AKS466",
  "DMSO_new" = "DMSO",
  "HPA-12_new" = "HPA-12",
  "Myriocin_new" = "Myriocin"
)
```

If needed, update the PCA plot colors here:

```r
scale_fill_manual(values = c(
  "AKS466" = "#5D3A9B",
  "DMSO" = "#E66100",
  "HPA-12" = "#0C7BDC",
  "Myriocin" = "#FFC20A"
))
```

## What the script does

### Data preparation and cleaning

* reads all Excel files in the input directory
* cleans column names
* assigns inhibitor identity from the file name
* aligns files with differing columns by taking the union of all column names
* combines all files into one data frame
* encodes `infection` and `inhibitor` as factors
* optionally subsets to uninfected samples only
* identifies lipid columns and coerces them to numeric
* applies `log10(x + 1)` transformation to lipid values

### Quality control filtering

* removes lipids with only missing values
* removes lipids with zero variance across samples
* removes samples with incomplete data
* z-scores lipids when Euclidean distance is used
* creates:

  * `df_cc` for filtered sample metadata
  * `X_cc` for the filtered lipid matrix

### PCA

* runs PCA on the filtered lipid matrix
* calculates variance explained by each principal component
* writes variance explained to a `.csv` file
* creates a PCA plot stored as `p_pca_2`

### PERMANOVA and dispersion checks

* runs PERMANOVA with `adonis2` to test the effect of `inhibitor`
* uses 999 permutations
* performs dispersion checks with `betadisper` and `permutest`
* writes the main PERMANOVA results and dispersion test results to text files

## Output

The script writes the following output files to the output directory:

* `pca_variance_explained.csv`
* `inhibitors_total_permanova_main.txt`
* `inhibitors_total_permanova_dispersion_tests.txt`

It also creates the following R objects in the session:

* `df`
* `df_cc`
* `X_cc`
* `pca_fit`
* `pca_var_df`
* `scores`
* `p_pca_2`
* `adon_uninfected_ind_dmso`
* `per_inhib`

## Notes

* The script assumes inhibitor identity can be assigned from workbook names.
* If workbook names differ from the current mapping, update `name_map`.
* Non-numeric lipid values may be converted to `NA`.
* Euclidean distance is appropriate here because the lipid data are z-scored before analysis.
* If `KEEP_ONLY_UNINFECTED <- TRUE`, infection-specific comparisons are not performed because only uninfected samples are retained.
* PCA, PERMANOVA, and dispersion checks are all performed within the same `.R` file.


