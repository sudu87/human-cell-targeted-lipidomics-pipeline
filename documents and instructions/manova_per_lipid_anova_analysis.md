# MANOVA + Per-Lipid ANOVA + Post-hoc FDR Workflow

This workflow tests whether the lipid profile differs across infection status, inhibitor treatment, and their interaction. It combines:

1. **MANOVA** as the global multivariate test across all lipid species.
2. **Per-lipid ANOVA** to identify which individual lipids contribute to group differences.
3. **Estimated marginal means post-hoc contrasts** for biologically targeted comparisons.
4. **Benjamini–Hochberg false discovery rate correction** for multiple testing.

The intended input is an Excel file containing one row per independent biological sample and one column per lipid species.

The script assumes two metadata columns:

| Column | Meaning |
|---|---|
| `infection` | Infection status, e.g. `Uninfected`, `infected` |
| `inhibitor` | Treatment group, e.g. `AKS466`, `HPA-12`, `Desipramin`, `Myriocin`, `ARC39`, `untreated`, `DMSO` |

All remaining numeric columns are treated as lipid abundance variables unless explicitly added to `additional_meta_cols` in the script.

---

Statistical model

The multivariate model is:

```r
data matrix ~ infection * inhibitor
```

This tests three terms:

| Term | Interpretation |
|---|---|
| `infection` | Overall lipidome difference between infected and uninfected samples |
| `inhibitor` | Overall lipidome difference among inhibitor/treatment groups |
| `infection:inhibitor` | Whether infection effects depend on inhibitor treatment |

The script reports MANOVA using **Pillai's trace** as the primary multivariate test statistic.

---

## Why Pillai's trace?

Pillai's trace is one of the standard MANOVA test statistics, alongside Wilks' lambda, Hotelling–Lawley trace, and Roy's largest root. The R `summary.manova()` method supports all four and uses Pillai as the first/default option. Pillai's trace is often preferred as a robust default when biological data may show moderate deviations from ideal MANOVA assumptions.

Interpretation:

- Higher Pillai's trace indicates stronger multivariate group separation.
- A small p-value, for example `p < 0.05`, indicates evidence against the null hypothesis that the multivariate lipid mean vectors are equal for that model term.
- A significant MANOVA result does **not** identify which individual lipid differs; follow-up per-lipid ANOVA and post-hoc contrasts are required.

---

## When MANOVA is appropriate

MANOVA is appropriate when:

- Multiple lipid species are measured from the same experimental unit.
- The lipid variables are biologically correlated.
- The goal is to test an overall lipidome-level shift.
- The samples are independent biological replicates.
- The number of response variables is reasonable relative to the residual degrees of freedom.

This is often conceptually appropriate for sphingolipid or lipidomics data because lipid species within the same pathway are typically correlated and are not independent biological endpoints.

---

## When MANOVA is not recommended

Avoid or interpret MANOVA cautiously when:

- The main goal is only to interpret each lipid individually.
- The number of lipids is larger than, or close to, the residual degrees of freedom.
- The design is very unbalanced or has missing infection-by-inhibitor cells.
- There are very few biological replicates per group.
- Many lipid variables are nearly perfectly correlated, causing rank-deficient residual covariance matrices.
- Distributions are highly skewed even after transformation.

For high-dimensional lipidomics data, consider reducing the response matrix before classical MANOVA, for example by filtering lipid species, using pathway-level summaries, or testing principal components. A permutation-based multivariate method such as PERMANOVA may also be appropriate for exploratory global lipidome testing, but that is a different null hypothesis and should be reported separately.

---

## Transformation

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

## Multiple testing and FDR strategy

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

## Post-hoc contrast families

The script performs two post-hoc families.

### A. Infection effect within each inhibitor

Question:

> Within each inhibitor, does infection change each lipid?

Model estimated marginal means:

```r
emmeans(fit, ~ infection | inhibitor)
```

### B. Inhibitor/treatment effects within each infection status

Question:

> Within infected or uninfected samples, which inhibitors differ from each other?

Model estimated marginal means:

```r
emmeans(fit, ~ inhibitor | infection)
```


