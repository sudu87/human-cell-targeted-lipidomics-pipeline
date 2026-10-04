library(readxl)
library(dplyr)
library(ggplot2)

# Standalone HeLa/HUVEC pairwise FDR heatmap.
# Configure paths directly below or override them with environment variables.
input_file <- Sys.getenv(
  "CELL_LINE_PAIRWISE_INPUT",
  unset = "path/to/hela_huvec_pairwise_results.xlsx"
)
output_dir <- Sys.getenv(
  "CELL_LINE_HEATMAP_OUTPUT_DIR",
  unset = "outputs/cell_line_pairwise_fdr_heatmap"
)
output_prefix <- Sys.getenv(
  "CELL_LINE_HEATMAP_PREFIX",
  unset = "hela_huvec_pairwise_fdr_heatmap"
)

alpha <- 0.05
fill_limits <- c(-1.5, 1.5)

total_lipid_order <- c(
  "dh_sph", "sph", "s1p",
  "dh_cer_total", "cer_total",
  "dh_sm_total", "sm_total",
  "hex_cer_total", "lac_cer_total"
)

individual_lipid_order <- c(
  "dh_sph", "sph", "s1p",
  "dh_cer_16_0", "dh_cer_18_0", "dh_cer_20_0",
  "dh_cer_22_0", "dh_cer_24_0", "dh_cer_24_1",
  "cer_16_0", "cer_18_0", "cer_20_0", "cer_22_0", "cer_24_0", "cer_24_1",
  "dh_sm_16_0", "dh_sm_18_0", "dh_sm_20_0",
  "dh_sm_22_0", "dh_sm_24_0", "dh_sm_24_1",
  "sm_16_0", "sm_18_0", "sm_20_0", "sm_22_0", "sm_24_0", "sm_24_1",
  "hex_cer_16_0", "hex_cer_24_1",
  "lac_cer_16_0", "lac_cer_24_1"
)

if (!file.exists(input_file)) {
  stop(
    "Input file does not exist. Set CELL_LINE_PAIRWISE_INPUT or update input_file: ",
    input_file
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df_heat <- read_excel(input_file)
required_cols <- c("lipid", "cell_line", "contrast", "estimate", "p_fdr")
missing_cols <- setdiff(required_cols, names(df_heat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_heat <- df_heat |>
  mutate(
    estimate = suppressWarnings(as.numeric(estimate)),
    p_fdr = suppressWarnings(as.numeric(p_fdr))
  ) |>
  filter(
    !is.na(lipid),
    !is.na(cell_line),
    !is.na(contrast),
    !is.na(estimate),
    !is.na(p_fdr)
  )

if (nrow(df_heat) == 0) {
  stop("No complete rows remain after input validation.")
}

duplicate_cells <- df_heat |>
  count(lipid, cell_line, contrast) |>
  filter(n > 1)
if (nrow(duplicate_cells) > 0) {
  stop("Duplicate lipid/cell-line/contrast combinations were found in the input.")
}

lipid_order <- if (any(grepl("_total$", df_heat$lipid))) {
  total_lipid_order
} else {
  individual_lipid_order
}

missing_from_order <- setdiff(unique(df_heat$lipid), lipid_order)
if (length(missing_from_order) > 0) {
  warning(
    "Lipids not listed in the configured order will be appended: ",
    paste(missing_from_order, collapse = ", ")
  )
}
lipid_levels <- c(
  lipid_order[lipid_order %in% unique(df_heat$lipid)],
  sort(missing_from_order)
)

plot_df <- df_heat |>
  mutate(
    log2FC = log2(10^estimate),
    significant = p_fdr < alpha,
    log2FC_sig = ifelse(significant, log2FC, NA_real_),
    outside_fill_limits = significant &
      (log2FC < fill_limits[1] | log2FC > fill_limits[2]),
    displayed_log2FC = ifelse(
      significant,
      pmax(fill_limits[1], pmin(fill_limits[2], log2FC)),
      NA_real_
    ),
    lipid = factor(lipid, levels = rev(lipid_levels)),
    cell_line = recode(
      tolower(as.character(cell_line)),
      "hela" = "HeLa",
      "huvec" = "HUVEC",
      .default = as.character(cell_line)
    ),
    cell_line = factor(cell_line, levels = c("HeLa", "HUVEC"))
  )

if (any(is.na(plot_df$cell_line))) {
  stop("Unexpected cell_line labels. Expected hela/HeLa and huvec/HUVEC.")
}

write.csv(
  plot_df,
  file.path(output_dir, paste0(output_prefix, "_plot_data.csv")),
  row.names = FALSE
)

audit_summary <- plot_df |>
  summarise(
    total_cells = n(),
    significant_cells = sum(significant),
    nonsignificant_cells = sum(!significant),
    significant_outside_fill_limits = sum(outside_fill_limits),
    fill_limit_min = fill_limits[1],
    fill_limit_max = fill_limits[2]
  )

write.csv(
  audit_summary,
  file.path(output_dir, paste0(output_prefix, "_scale_audit.csv")),
  row.names = FALSE
)

p_heatmap <- ggplot(
  plot_df,
  aes(x = contrast, y = lipid, fill = log2FC_sig)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_y_discrete(position = "right", drop = FALSE) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = fill_limits,
    breaks = seq(fill_limits[1], fill_limits[2], by = 0.5),
    oob = scales::squish,
    na.value = "grey85",
    name = "log2 fold-change"
  ) +
  facet_grid(. ~ cell_line) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 0),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.title = element_text(size = 9)
  )

ggsave(
  filename = file.path(output_dir, paste0(output_prefix, ".pdf")),
  plot = p_heatmap,
  width = 4.2,
  height = if (length(lipid_levels) > 15) 7 else 4.5,
  units = "in"
)

message(
  "Saved heatmap. Significant cells outside the fill limits were saturated, not greyed: ",
  audit_summary$significant_outside_fill_limits
)
