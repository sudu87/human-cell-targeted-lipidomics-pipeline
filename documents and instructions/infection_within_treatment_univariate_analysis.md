# Compare infected vs uninfected within one treatment

This document explains how to use the script `scripts/infection_within_treatment_univariate_analysis.R`.

## Instructions

- Place your input Excel file in your project folder or note its full file path.
- Open `scripts/infection_within_treatment_univariate_analysis.R`.
- Replace the placeholder path in the following line with the path to your own Excel file:

```r
file_path <- "path/to/your/input_file.xlsx"
````

* Replace the placeholder output path in the following line with the file path where you want the results to be saved:

```r
output_file <- "path/to/your/output_file.xlsx"
```

* Make sure your Excel file is an `.xlsx` file.
* Make sure your dataset contains a column named `infection`.
* Arrange all lipid species measurements in the remaining columns.
* Check that the `infection` column contains the exact labels `Uninfected` and `Infected`.
* Keep samples as separate rows.
* Make sure lipid measurement columns contain numeric values, because non-numeric entries may be converted to `NA`.
* The script applies a `log10(x + 1)` transformation before ANOVA.
* Run the script after loading the required R packages.
* Use the resulting object `res_infection` for later inspection or export.

## Overview

This script compares infected and uninfected samples within one treatment condition using one-way ANOVA for each lipid. Pairwise comparisons are then calculated using `emmeans`, and p-values are adjusted using FDR correction.

## Input requirements

Your Excel file must:

* be provided as an `.xlsx` file
* contain a metadata column named `infection`
* contain lipid species measurements in all remaining columns
* keep samples as separate rows
* use the expected infection labels

Expected values:

* `infection`: `Uninfected`, `Infected`

## Specific lines you may need to edit

Update this line with the path to your own Excel file:

```r
file_path <- "path/to/your/input_file.xlsx"
```

Update this line with the path to your output file:

```r
output_file <- "path/to/your/output_file.xlsx"
```

If needed, edit the infection factor order here:

```r
infection = factor(infection, levels = c("Uninfected", "Infected"))
```

## What the script does

The script:

* reads the Excel file
* standardizes column names
* cleans and formats the `infection` column
* identifies lipid columns
* converts lipid columns to numeric values
* applies `log10(x + 1)` transformation before analysis
* runs one-way ANOVA for each lipid using `infection` as the predictor
* calculates pairwise comparisons with `emmeans`
* adjusts p-values using FDR correction
* combines all lipid-wise results into the object `res_infection`
* writes the results to an Excel file

## Output

The final statistical results are stored in:

```r
res_infection
```

The results are also written to the Excel file defined in:

```r
output_file
```

## Notes

* Non-numeric values in lipid columns may be converted to `NA`.
* Zero values are allowed because the script uses `log10(x + 1)`.
* The script assumes that only one treatment condition is present in the input file.
* Make sure the required packages are loaded before running the script.



