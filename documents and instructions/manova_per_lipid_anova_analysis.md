# Total Lipids Dataset: MANOVA + Per-Lipid ANOVA + Post-hoc FDR Workflow

This workflow tests whether the **total lipid profile** differs across infection status, inhibitor treatment, and their interaction. It combines:

1. **MANOVA** as the global multivariate test across all lipid species.
2. **Per-lipid ANOVA** to identify which individual lipids contribute to group differences.
3. **Estimated marginal means post-hoc contrasts** for biologically targeted comparisons.
4. **Benjamini–Hochberg false discovery rate correction** for multiple testing.

The intended input is an Excel file containing one row per independent biological sample and one column per lipid species.

---

## 1. Scientific question

Use this analysis when the main biological question is:

> Does infection, inhibitor treatment, or their interaction shift the lipidome as a multivariate profile?

For example, with sphingolipid species measured from infected and uninfected cells across several inhibitors, MANOVA tests whether the combined lipid abundance vector differs among experimental groups.

---

## 2. Example experimental design

The script assumes two metadata columns:

| Column | Meaning |
|---|---|
| `infection` | Infection status, e.g. `Uninfected`, `infected` |
| `inhibitor` | Treatment group, e.g. `AKS466`, `HPA-12`, `Desipramin`, `Myriocin`, `ARC39`, `untreated`, `DMSO` |

All remaining numeric columns are treated as lipid abundance variables unless explicitly added to `additional_meta_cols` in the script.

---

## 3. Statistical model

The multivariate model is:

```r
cbind(lipid_1, lipid_2, ..., lipid_p) ~ infection * inhibitor
```

This tests three terms:

| Term | Interpretation |
|---|---|
| `infection` | Overall lipidome difference between infected and uninfected samples |
| `inhibitor` | Overall lipidome difference among inhibitor/treatment groups |
| `infection:inhibitor` | Whether infection effects depend on inhibitor treatment |

The script reports MANOVA using **Pillai's trace** as the primary multivariate test statistic.

---

## 4. Why Pillai's trace?

Pillai's trace is one of the standard MANOVA test statistics, alongside Wilks' lambda, Hotelling–Lawley trace, and Roy's largest root. The R `summary.manova()` method supports all four and uses Pillai as the first/default option. Pillai's trace is often preferred as a robust default when biological data may show moderate deviations from ideal MANOVA assumptions.

Interpretation:

- Higher Pillai's trace indicates stronger multivariate group separation.
- A small p-value, for example `p < 0.05`, indicates evidence against the null hypothesis that the multivariate lipid mean vectors are equal for that model term.
- A significant MANOVA result does **not** identify which individual lipid differs; follow-up per-lipid ANOVA and post-hoc contrasts are required.

---

## 5. When MANOVA is appropriate

MANOVA is appropriate when:

- Multiple lipid species are measured from the same experimental unit.
- The lipid variables are biologically correlated.
- The goal is to test an overall lipidome-level shift.
- The samples are independent biological replicates.
- The number of response variables is reasonable relative to the residual degrees of freedom.

This is often conceptually appropriate for sphingolipid or lipidomics data because lipid species within the same pathway are typically correlated and are not independent biological endpoints.

---

## 6. When MANOVA is not recommended

Avoid or interpret MANOVA cautiously when:

- The main goal is only to interpret each lipid individually.
- The number of lipids is larger than, or close to, the residual degrees of freedom.
- The design is very unbalanced or has missing infection-by-inhibitor cells.
- There are very few biological replicates per group.
- Many lipid variables are nearly perfectly correlated, causing rank-deficient residual covariance matrices.
- Distributions are highly skewed even after transformation.

For high-dimensional lipidomics data, consider reducing the response matrix before classical MANOVA, for example by filtering lipid species, using pathway-level summaries, or testing principal components. A permutation-based multivariate method such as PERMANOVA may also be appropriate for exploratory global lipidome testing, but that is a different null hypothesis and should be reported separately.

---

## 7. Transformation

The script applies:

```r
log10(x + 1)
```

by default.

This is useful for lipid abundance data because it reduces right-skew and compresses very large abundance differences. If the values are already log-transformed, set:

```r
log10_transform <- FALSE
```

in the script.

---

## 8. Multiple testing and FDR strategy

The post-hoc contrasts are first estimated using `emmeans` with unadjusted pairwise tests. After all lipids and contrasts are combined, the script applies Benjamini–Hochberg FDR correction.

The script reports:

| Column | Meaning |
|---|---|
| `p.value` | Raw p-value from the contrast |
| `p_adj_bh_by_lipid` | BH correction within each lipid for that contrast family |
| `p_adj_bh_global` | BH correction across all lipids and contrasts in that post-hoc family |
| `significant_global_fdr_0_05` | Whether `p_adj_bh_global < 0.05` |

For publication, use `p_adj_bh_global` unless there is a pre-specified reason to control FDR separately within each lipid.

---

## 9. Post-hoc contrast families

The script performs two post-hoc families.

### A. Infection effect within each inhibitor

Question:

> Within each inhibitor, does infection change each lipid?

Model estimated marginal means:

```r
emmeans(fit, ~ infection | inhibitor)
```

Output file:

```text
myriocin_posthoc_infection_contrasts_by_inhibitor_bh_fdr.xlsx
```

### B. Inhibitor/treatment effects within each infection status

Question:

> Within infected or uninfected samples, which inhibitors differ from each other?

Model estimated marginal means:

```r
emmeans(fit, ~ inhibitor | infection)
```

Output file:

```text
myriocin_posthoc_inhibitor_contrasts_by_infection_bh_fdr.xlsx
```

---

## 10. Output files

The script creates an output folder:

```text
statistics/total_lipid_stats/myriocin/
```

Main outputs:

| File | Contents |
|---|---|
| `myriocin_statistics_results.xlsx` | Main workbook containing MANOVA, ANOVA, post-hoc tests, group counts, residual diagnostics |
| `myriocin_manova_pillai_summary.txt` | Plain-text MANOVA Pillai summary |
| `myriocin_univariate_aov_summary.txt` | Plain-text per-lipid ANOVA summary from the MANOVA fit |
| `myriocin_posthoc_infection_contrasts_by_inhibitor_bh_fdr.xlsx` | Infection contrasts within each inhibitor |
| `myriocin_posthoc_inhibitor_contrasts_by_infection_bh_fdr.xlsx` | Inhibitor contrasts within each infection group |
| `myriocin_processed_data_log10.xlsx` | Processed data used for modelling |
| `myriocin_session_info.txt` | R session information for reproducibility |
| `plots/pooled_standardized_residual_density.png` | Residual density diagnostic |
| `plots/pooled_standardized_residual_qq.png` | Residual Q-Q diagnostic |

If the optional `biotools` package is installed, the script also writes a Box's M covariance homogeneity diagnostic.

---

## 11. How to run

Recommended GitHub structure:

```text
project/
├── averaged_sphingolipidspecies/
│   └── Myriocin.xlsx
├── scripts/
│   └── total_lipids_manova_analysis.R
└── statistics/
    └── total_lipid_stats/
```

Run from the project root:

```bash
Rscript scripts/total_lipids_manova_analysis.R
```

Install required R packages:

```r
install.packages(c(
  "readxl",
  "janitor",
  "dplyr",
  "tidyr",
  "purrr",
  "broom",
  "emmeans",
  "writexl",
  "ggplot2",
  "tibble"
))
```

Optional diagnostic package:

```r
install.packages("biotools")
```

---

## 12. Reporting template for manuscript

Example wording:

> Lipid abundance values were log10(x + 1)-transformed before statistical analysis. Global lipidome differences were tested using MANOVA with infection, inhibitor treatment, and their interaction as fixed factors. Pillai's trace was used as the primary multivariate test statistic. For individual lipid species, follow-up two-way ANOVAs were fitted using the same model structure. Post-hoc comparisons were estimated using estimated marginal means. P-values were adjusted for multiple testing using the Benjamini–Hochberg false discovery rate procedure.

When reporting results, include:

- Number of biological replicates per infection-by-inhibitor group.
- Number of lipid variables included in MANOVA.
- Transformation used.
- MANOVA Pillai's trace, approximate F statistic, numerator df, denominator df, and p-value.
- FDR correction method and family of tests corrected.
- Whether post-hoc estimates are reported on transformed scale or back-transformed ratio scale.

---

## 13. Important interpretation notes

- MANOVA tests a **global multivariate null hypothesis**; it does not identify specific lipids.
- Per-lipid ANOVA and post-hoc tests are follow-up analyses and require multiple-testing correction.
- If the interaction term `infection:inhibitor` is significant, interpret simple effects first, such as infection within each inhibitor or inhibitor contrasts within each infection status.
- The post-hoc `estimate` column is on the transformed scale when `log10_transform <- TRUE`.
- The optional `ratio_log10p1` column is `10^estimate`, which approximates a ratio on the `x + 1` scale, not a raw untransformed fold-change.

---

## 14. Sources

- R Core Team. `summary.manova`: Summary Method for Multivariate Analysis of Variance. https://stat.ethz.ch/R-manual/R-devel/library/stats/html/summary.manova.html
- Penn State STAT 505. Lesson 8: Multivariate Analysis of Variance (MANOVA). https://online.stat.psu.edu/stat505/book/export/html/762
- Lenth RV. `emmeans` package documentation and quick-start guide. https://rvlenth.github.io/emmeans/articles/AQuickStart.html
- R Core Team. `p.adjust`: Adjust P-values for Multiple Comparisons. https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
- Benjamini Y, Hochberg Y. 1995. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society Series B 57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Olson CL. 1974. Comparative Robustness of Six Tests in Multivariate Analysis of Variance. Journal of the American Statistical Association 69:894–908. https://doi.org/10.1080/01621459.1974.10480224
