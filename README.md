# Targeted lipidomics analysis of sphingolipid metabolism during _Simkania negevensis_ infection

![Lipidomics workflow overview](images/lipidomics3_250526.png)

## Related publication

**Chlamydia-like bacterium Simkania negevensis exploits host sphingolipid salvage pathway and sphingomyelin synthesis during infection**

Mohanty, A., Weinrich, J. D., Schumacher, F., Rühling, M., Sunuwar, S., Szegedi, H., Wigger, D., Schmelz, F., Panda, B. K., Kappe, C., Brenner, D., Schirmer, M., Arenz, C., Seibel, J., Holthuis, J. C. M., Das, S., Fraunholz, M., Kleuser, B., and Kozjak-Pavlovic, V.

## Overview

This repository contains R scripts and documentation for the analysis of targeted lipidomics data from human cells infected with _Simkania negevensis_.

The analysis workflow includes:

- data import and preprocessing
- data cleaning and quality assessment
- replicate averaging
- lipid abundance visualization
- heatmap generation
- multivariate statistical analysis
- PERMANOVA and dispersion testing
- downstream exploratory visualization

## Repository structure

```text
├─ README.md
├─ REPRODUCIBILITY.md
├─ CITATION.cff
├─ LICENSE
├─ demo_data/
│  ├─ README.md
│  └─ sms12_demo_lipidomics.xlsx
├─ documents and instructions/
│  ├─ transformation_diagnostics.md
│  ├─ sms12_manova_pca_dispersion_analysis.md
│  ├─ sms12_individual_lipids_permanova_dispersion_analysis.md
│  ├─ custom_pairwise_lipid_contrast_heatmap.md
│  ├─ pairwise_fdr_lipid_heatmap.md
│  └─ ...
├─ scripts/
│  ├─ transformation_diagnostics.R
│  ├─ sms12_manova_pca_dispersion_analysis.R
│  ├─ sms12_individual_lipids_permanova_dispersion_analysis.R
│  ├─ custom_pairwise_lipid_contrast_heatmap.R
│  ├─ pairwise_fdr_lipid_heatmap.R
│  └─ ...
└─ images/
   ├─ hist_raw_values.png
   ├─ hist_log10_values.png
   ├─ density_raw_vs_log.png
   └─ qqplot_residuals.png
```

## Getting started

Before running the analysis scripts, make sure that R and the required packages are installed.

Required CRAN packages used across the scripts include:

```r
readxl
janitor
dplyr
stringr
purrr
tidyr
tibble
ggplot2
pheatmap
vegan
emmeans
broom
writexl
rcompanion
scales
```

For a reproducible package environment, use the versions recorded in [`renv.lock`](renv.lock):

```r
install.packages("renv")
renv::restore()
```

Alternatively, install missing packages manually with:

```r
install.packages(c(
  "readxl",
  "janitor",
  "dplyr",
  "stringr",
  "purrr",
  "tidyr",
  "tibble",
  "ggplot2",
  "pheatmap",
  "vegan",
  "emmeans",
  "broom",
  "writexl",
  "rcompanion",
  "scales"
))
```

The `tools` package ships with R and is loaded by scripts that need it.

## Running the analysis

Scripts are located in the `scripts/` directory. Each script performs a specific part of the lipidomics workflow.

Each script is intended to run as a standalone analysis file after you install the required packages listed above and update the user-defined input/output paths near the top of the script. The scripts load their own libraries with `library()` calls, so you do not need to source another project file first.

For example:

```r
source("scripts/transformation_diagnostics.R")
source("scripts/pairwise_fdr_lipid_heatmap.R")
source("scripts/sms12_manova_pca_dispersion_analysis.R")
```

Detailed explanations of selected workflows are provided in the `documents and instructions/` directory.

For step-by-step reproducibility instructions, including Zenodo input-file mapping and example run commands, see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

A small synthetic demo workbook is provided in [`demo_data/`](demo_data/) for installation checks and reviewer testing.

## Input data availability

Primary experimental input data files are not included in this repository. They are deposited in Zenodo:

Das, S., & Mohanty, A. (2026). *Targeted lipidomics analysis of sphingolipid metabolism during Simkania negevensis infection* (v.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.18866967

The Zenodo record is currently under embargo. Access to files may be restricted until the embargo is lifted.

Please ensure that file paths inside the scripts are adjusted to match your local data directory.

## Outputs

The scripts generate exploratory plots, heatmaps, and statistical summaries for lipid abundance patterns across experimental conditions.

Typical outputs include:

- raw and log-transformed value distributions
- replicate-averaged heatmaps
- pairwise comparison and fdr heatmaps
- PERMANOVA results
- dispersion test results
- diagnostic plots

## Contact

For questions about the lipidomics analysis pipeline, scripts, or reproducibility, please contact:

**Sudip Das**  
Email: sudip.das@tum.de

**Arpita Mohanty**  
Email: arpita.mohanty@uni-wuerzburg.de

## Citation

If you use this repository or adapt the analysis workflow, please cite the related publication. Citation metadata is available in [`CITATION.cff`](CITATION.cff), and GitHub will show a “Cite this repository” option for it.

## License

This repository is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). See [`LICENSE`](LICENSE) for details.
