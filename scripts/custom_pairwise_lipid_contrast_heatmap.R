library(readxl)
library(dplyr)
library(ggplot2)

# ==============================
# Custom pairwise lipid contrast heatmap
# - Reads a pairwise lipid contrast results table
# - Keeps user-selected contrasts in a user-defined x-axis order
# - Keeps total lipid classes or individual lipid species in a user-defined y-axis order
# - Plots significant log2 fold-changes after FDR filtering
# ==============================

## ---- Configure user-defined input and output paths ----
input_file <- "path/to/your/pairwise_results.xlsx"
sheet_name <- NULL
output_dir <- "path/to/your/output_directory"
output_prefix <- "custom_pairwise_lipid_contrast_heatmap"

## ---- Configure input columns ----
lipid_col <- "lipid"
contrast_col <- "contrast"
estimate_col <- "estimate"
p_value_col <- "p_fdr"
infection_col <- "infection"

## ---- Configure contrast and lipid ordering ----
# This vector chooses which pairwise comparisons are shown and sets the x-axis order.
keep_contrasts <- c(
  "d_sms12_DKO - h_sms12_WT",
  "a_sms1_AHT - e_sms1_no_AHT",
  "e_sms1_no_AHT - h_sms12_WT",
  "a_sms1_AHT - h_sms12_WT",
  "b_sms2_AHT - f_sms2_no_AHT",
  "f_sms2_no_AHT - h_sms12_WT",
  "b_sms2_AHT - h_sms12_WT",
  "(c_sms2-M64R_AHT) - (g_sms2-M64R_no_AHT)",
  "(g_sms2-M64R_no_AHT) - h_sms12_WT",
  "(c_sms2-M64R_AHT) - h_sms12_WT"
)

# This vector chooses which lipids are shown and sets the y-axis order.
# It can contain total lipid classes or individual lipid species.
custom_lipid_order <- c(
  "dh_sph", "sph", "s1p",
  "dh_cer16_0", "dh_cer18_0", "dh_cer20_0", "dh_cer22_0", "dh_cer24_0", "dh_cer24_1",
  "cer16_0", "cer18_0", "cer20_0", "cer22_0", "cer24_0", "cer24_1",
  "dh_sm16_0", "dh_sm18_0", "dh_sm20_0", "dh_sm22_0", "dh_sm24_0", "dh_sm24_1",
  "sm16_0", "sm18_0", "sm20_0", "sm22_0", "sm24_0", "sm24_1",
  "hex_cer16_0", "hex_cer24_1",
  "lac_cer16_0", "lac_cer24_1"
)

keep_only_custom_lipids <- TRUE

## ---- Configure plotting options ----
alpha <- 0.05
estimate_scale <- "log10"
fill_limits <- c(-8, 8)
infection_order <- c("no", "yes")
infection_labels <- c("no" = "Uninfected", "yes" = "Infected")
show_x_axis_labels <- FALSE

## ---- Read input ----
if (!file.exists(input_file)) {
  stop("Input file does not exist. Update input_file: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df_heat <- if (is.null(sheet_name) || identical(sheet_name, "")) {
  read_excel(input_file)
} else {
  read_excel(input_file, sheet = sheet_name)
}

required_cols <- c(lipid_col, contrast_col, estimate_col, p_value_col, infection_col)
missing_cols <- setdiff(required_cols, names(df_heat))
if (length(missing_cols) > 0) {
  stop("Missing required input columns: ", paste(missing_cols, collapse = ", "))
}

df_heat <- df_heat |>
  rename(
    lipid = all_of(lipid_col),
    contrast = all_of(contrast_col),
    estimate = all_of(estimate_col),
    p_value = all_of(p_value_col),
    infection = all_of(infection_col)
  ) |>
  mutate(
    estimate = suppressWarnings(as.numeric(estimate)),
    p_value = suppressWarnings(as.numeric(p_value))
  )

## ---- Helpers ----
estimate_to_log2fc <- function(x, scale) {
  if (scale == "log10") {
    return(log2(10^x))
  }

  if (scale == "log2") {
    return(x)
  }

  if (scale == "natural_log") {
    return(x / log(2))
  }

  if (scale == "fold_change") {
    return(log2(x))
  }

  stop(
    "Unknown estimate_scale: ",
    scale,
    ". Use one of: log10, log2, natural_log, fold_change."
  )
}

## ---- Filter and transform ----
missing_contrasts <- setdiff(keep_contrasts, unique(df_heat$contrast))
if (length(missing_contrasts) > 0) {
  warning(
    "These keep_contrasts were not found in the input and will be absent: ",
    paste(missing_contrasts, collapse = ", ")
  )
}

missing_lipids <- setdiff(custom_lipid_order, unique(df_heat$lipid))
if (length(missing_lipids) > 0) {
  warning(
    "These custom_lipid_order entries were not found in the input: ",
    paste(missing_lipids, collapse = ", ")
  )
}

df_hm <- df_heat |>
  filter(contrast %in% keep_contrasts) |>
  filter(
    !is.na(lipid),
    !is.na(estimate),
    !is.na(p_value),
    !is.na(infection),
    !is.na(contrast)
  )

if (keep_only_custom_lipids) {
  df_hm <- df_hm |>
    filter(lipid %in% custom_lipid_order)
}

if (nrow(df_hm) == 0) {
  stop("No rows remain after filtering. Check keep_contrasts and custom_lipid_order.")
}

custom_present <- custom_lipid_order[custom_lipid_order %in% unique(df_hm$lipid)]
other_lipids <- setdiff(sort(unique(df_hm$lipid)), custom_present)
lipid_levels <- if (keep_only_custom_lipids) {
  custom_present
} else {
  c(custom_present, other_lipids)
}

contrast_levels <- keep_contrasts[keep_contrasts %in% unique(df_hm$contrast)]

extra_infections <- setdiff(unique(as.character(df_hm$infection)), infection_order)
if (length(extra_infections) > 0) {
  warning(
    "Found infection values not listed in infection_order; appending them: ",
    paste(extra_infections, collapse = ", ")
  )
  infection_order <- c(infection_order, extra_infections)
  infection_labels <- c(infection_labels, setNames(extra_infections, extra_infections))
}

plot_df <- df_hm |>
  mutate(
    log2FC = estimate_to_log2fc(estimate, estimate_scale),
    log2FC_sig = ifelse(p_value < alpha, log2FC, NA_real_),
    lipid = factor(lipid, levels = rev(lipid_levels)),
    contrast = factor(contrast, levels = contrast_levels),
    infection = recode(as.character(infection), !!!infection_labels),
    infection = factor(infection, levels = unname(infection_labels[infection_order]))
  )

write.csv(
  plot_df,
  file = file.path(output_dir, paste0(output_prefix, "_plot_data.csv")),
  row.names = FALSE
)

writeLines(
  missing_contrasts,
  con = file.path(output_dir, paste0(output_prefix, "_missing_contrasts.txt"))
)

writeLines(
  missing_lipids,
  con = file.path(output_dir, paste0(output_prefix, "_missing_lipids.txt"))
)

## ---- Heatmap ----
x_axis_text <- if (show_x_axis_labels) {
  element_text(angle = 45, hjust = 1, size = 8)
} else {
  element_blank()
}

p_heatmap <- ggplot(plot_df, aes(x = contrast, y = lipid, fill = log2FC_sig)) +
  scale_y_discrete(position = "right", drop = FALSE) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = fill_limits,
    na.value = "grey85",
    guide = guide_colorbar(barwidth = 0.7)
  ) +
  facet_grid(. ~ infection, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 0),
    axis.text.x = x_axis_text,
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(output_dir, paste0(output_prefix, ".pdf")),
  plot = p_heatmap,
  width = 10,
  height = 7,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, paste0(output_prefix, ".png")),
  plot = p_heatmap,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)
