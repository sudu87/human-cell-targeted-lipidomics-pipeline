
# PERMANOVA inhibitor analysis

This document explains how to use the script `scripts/permanova_inhibitor_analysis.R`.

## Instructions

- Place your input Excel files in one folder.
- Open `scripts/permanova_inhibitor_analysis.R`.
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
* Check that inhibitor identity can be assigned correctly from the file names.
* Update the file name mapping in the script if your workbook names differ from the current examples.
* Make sure the `infection` labels match the expected values used in the script.
* Make sure lipid measurement columns contain numeric values, because non-numeric values may be converted to `NA`.
* Run the script after loading the required R packages.
* Review the output text files written to the output directory.

## Overview

This script performs data preparation, quality control filtering, and PERMANOVA-based multivariate analysis for lipidomics datasets stored across multiple Excel files. Each workbook is treated as one inhibitor condition, the inhibitor label is assigned from the file name, and all files are combined into a single analysis table.

The script then:

* cleans and aligns the imported datasets
* encodes `infection` and `inhibitor` as categorical factors
* coerces lipid columns to numeric
* applies `log10(x + 1)` transformation
* removes lipids with only missing values
* removes lipids with zero variance
* removes samples with incomplete data
* z-scores lipid values when Euclidean distance is used
* runs PERMANOVA for global group effects
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

Adjust the analysis settings here if needed:

```r
METHOD <- "euclidean"
MIN_N  <- 4
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

If needed, edit the factor order here:

```r
df$infection <- factor(df$infection, levels = c("Uninfected", "Infected"))
```

## What the script does

### Data preparation and cleaning

* reads all Excel files in the input directory
* cleans column names
* assigns inhibitor identity from the file name
* aligns sheets with differing columns by taking the union of all column names
* combines all files into one data frame
* encodes `infection` and `inhibitor` as factors
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

### PERMANOVA and dispersion checks

* runs PERMANOVA with `adonis2` to test the effects of `infection` and `inhibitor`
* uses 999 permutations
* uses Euclidean distance on z-scored data when `METHOD = "euclidean"`
* performs dispersion checks with `betadisper` and `permutest`
* writes the main PERMANOVA results and dispersion test results to text files

## Output

The script writes the following output files to the output directory:

* `inhibitors_ind_permanova_main.txt`
* `inhibitors_ind_permanova_dispersion_tests.txt`

It also creates the following R objects in the session:

* `df`
* `df_cc`
* `X_cc`
* `adon_main`
* `per_inhib`
* `per_infec`

## Notes

* The script assumes inhibitor identity can be assigned from workbook names.
* If workbook names differ from the current mapping, update `name_map`.
* Non-numeric lipid values may be converted to `NA`.
* Euclidean distance is appropriate here because the lipid data are z-scored before analysis.
* The script keeps both data preparation and PERMANOVA steps in the same `.R` file.



