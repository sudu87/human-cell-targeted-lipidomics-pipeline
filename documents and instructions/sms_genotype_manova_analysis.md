# SMS genotype MANOVA analysis

This document explains how to use the script `scripts/sms_genotype_manova_analysis.R`.

## Instructions

- Place your input Excel file in your project folder.
- Open `scripts/sms_genotype_manova_analysis.R`.
- Replace the placeholder input file path with the path to your Excel file:

```r
input_file <- "path/to/your/input_file.xlsx"
```

- Replace the placeholder sheet name with the worksheet you want to analyze:

```r
sheet_name <- "your_sheet_name"
```

- Replace the placeholder output directory with the folder where you want the MANOVA results to be saved:

```r
out_dir <- "path/to/your/output_directory/"
```

- Make sure the Excel file contains metadata columns named `infection` and `condition`.
- Keep rows as individual samples and columns as metadata plus lipid measurements.
- Make sure lipid measurement columns contain numeric values.
- Run the script after loading the required R packages.
- Review the MANOVA output written to the output directory.

## Overview

This script performs MANOVA on lipidomics data from genotype or condition comparisons.

The script:

- reads one Excel worksheet
- cleans column names
- encodes `infection` and `condition` as categorical factors
- identifies lipid columns
- applies `log10(x + 1)` transformation to lipid values
- creates a lipid response matrix
- runs MANOVA using the model `infection * condition`
- reports Pillai's trace as the main MANOVA test

## Input requirements

Your Excel worksheet must:

- contain a metadata column named `infection`
- contain a metadata column named `condition`
- contain lipid measurements in the remaining columns
- have rows corresponding to individual samples
- have columns corresponding to metadata and lipid measurements

## Specific lines you may need to edit

Update the input file path:

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

Update the infection factor levels if your dataset uses different labels:

```r
infection = factor(infection, levels = c("no", "yes"))
```

## What the script does

### Data preparation

- reads the selected Excel worksheet
- cleans column names
- converts `infection` and `condition` to factors
- defines metadata columns
- treats all remaining columns as lipid measurements

### Transformation

- applies `log10(x + 1)` transformation to lipid measurements

### MANOVA

- creates a lipid response matrix
- runs MANOVA with the formula:

```r
resp_mat ~ infection * condition
```

- prints the MANOVA result using Pillai's trace
- saves the MANOVA result to the output directory

## Output

The script writes the following output file:

- `sms_genotype_manova_pillai.txt`

It also creates the following R objects in the session:

- `sms_data`
- `meta_cols`
- `lipid_cols`
- `resp_mat`
- `sms_manova_fit`
- `sms_manova_pillai`

## Notes

- `condition` is treated as the genotype or experimental condition variable.
- Pillai's trace is used as the main MANOVA statistic.
- The script assumes all non-metadata columns are lipid measurements.
- If your infection labels are not `no` and `yes`, update the factor levels in the script.
