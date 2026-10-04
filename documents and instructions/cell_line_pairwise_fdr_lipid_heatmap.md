# HeLa/HUVEC pairwise FDR lipid heatmap

This standalone workflow creates the HeLa/HUVEC infection-effect heatmap from a
pairwise results workbook. It supports either total lipid classes or individual
lipid species and preserves the configured biological lipid order.

## Input

The Excel workbook must contain:

- `lipid`
- `cell_line`
- `contrast`
- `estimate`
- `p_fdr`

`estimate` is assumed to come from a model fitted on a base-10 log scale. The
script converts it to log2 fold-change with `log2(10^estimate)`.

## Run

```sh
CELL_LINE_PAIRWISE_INPUT="path/to/hela_huvec_pairwise_results.xlsx" \
CELL_LINE_HEATMAP_OUTPUT_DIR="outputs/cell_line_pairwise_fdr_heatmap" \
Rscript scripts/cell_line_pairwise_fdr_lipid_heatmap.R
```

The script automatically selects the total-class or individual-species lipid
order based on the input lipid names.

## Significance and color limits

Cells with `p_fdr < 0.05` are colored by log2 fold-change. Non-significant cells
are grey. The color scale is fixed at -1.5 to +1.5 so that all effects in the
audited HeLa/HUVEC total and individual datasets are displayed without clipping.

Significant effects beyond these limits are saturated at the nearest endpoint
using `scales::squish`. This is important: ggplot2's default behavior censors
out-of-range values to `NA`, which would otherwise make a significant tile look
grey and therefore falsely non-significant.

## Outputs

- `hela_huvec_pairwise_fdr_heatmap.pdf`
- `hela_huvec_pairwise_fdr_heatmap_plot_data.csv`
- `hela_huvec_pairwise_fdr_heatmap_scale_audit.csv`

The plot-data file includes `outside_fill_limits` and `displayed_log2FC`. These
fields distinguish statistical significance from visual saturation at the
legend endpoints.
