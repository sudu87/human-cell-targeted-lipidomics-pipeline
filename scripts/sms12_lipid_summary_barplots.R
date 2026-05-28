library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================
# SMS1/2 lipid summary bar plots
# - Reads one raw Excel sheet containing infection, condition, and lipid columns
# - Summarises selected lipids by condition and infection
# - Plots mean +/- SE grouped bar plots for each selected lipid
# ==============================

## ---- Configure user-defined input and output paths ----
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"

## ---- Configure plot groups ----
# Set to NULL to use the order found in the input file.
condition_order <- c(
  "sms12_WT",
  "sms12_DKO",
  "sms1_16_no_AHT",
  "sms1_16_AHT",
  "sms2_no_AHT",
  "sms2_AHT",
  "sms2-M64R_no_AHT",
  "sms2-M64R_AHT"
)

infection_order <- c("no", "yes")
infection_labels <- c("no" = "Uninfected", "yes" = "Infected")

lipids_interest <- c(
  "cer_total",
  "sm_total",
  "hex_cer_total",
  "lac_cer_total",
  "dh_sm_total",
  "sph"
)

lipid_labels <- c(
  "cer_total" = "cer_total (pmol / sample)",
  "sm_total" = "SM (pmol / sample)",
  "hex_cer_total" = "hex_cer_total (pmol / sample)",
  "lac_cer_total" = "lac_cer_total (pmol / sample)",
  "dh_sm_total" = "dh_sm_total (pmol / sample)",
  "sph" = "sph_total (pmol / sample)"
)

if (!file.exists(input_file)) {
  stop("Input file does not exist. Update input_file: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Read and prepare raw data ----
df_sms12_raw <- read_excel(input_file, sheet = sheet_name) |>
  clean_names()

required_cols <- c("infection", "condition")
missing_required <- setdiff(required_cols, names(df_sms12_raw))
if (length(missing_required) > 0) {
  stop(
    "Missing required columns after clean_names(): ",
    paste(missing_required, collapse = ", ")
  )
}

missing_lipids <- setdiff(lipids_interest, names(df_sms12_raw))
if (length(missing_lipids) > 0) {
  warning(
    "These lipids were not found and will be skipped: ",
    paste(missing_lipids, collapse = ", ")
  )
}

lipids_present <- intersect(lipids_interest, names(df_sms12_raw))
if (length(lipids_present) == 0) {
  stop("None of the lipids in lipids_interest were found in the input file.")
}

extra_infections <- setdiff(
  unique(as.character(df_sms12_raw$infection)),
  infection_order
)
if (length(extra_infections) > 0) {
  warning(
    "Found infection values not listed in infection_order; appending them: ",
    paste(extra_infections, collapse = ", ")
  )
  infection_order <- c(infection_order, extra_infections)
}

if (is.null(condition_order)) {
  condition_order <- unique(as.character(df_sms12_raw$condition))
} else {
  extra_conditions <- setdiff(
    unique(as.character(df_sms12_raw$condition)),
    condition_order
  )

  if (length(extra_conditions) > 0) {
    warning(
      "Found condition values not listed in condition_order; appending them: ",
      paste(extra_conditions, collapse = ", ")
    )
    condition_order <- c(condition_order, extra_conditions)
  }
}

df_sms12_raw <- df_sms12_raw |>
  mutate(
    infection = factor(infection, levels = infection_order),
    condition = factor(condition, levels = condition_order)
  ) |>
  mutate(across(all_of(lipids_present), ~ suppressWarnings(as.numeric(.x))))

## ---- Reshape and summarise selected lipids ----
df_long <- df_sms12_raw |>
  select(infection, condition, all_of(lipids_present)) |>
  pivot_longer(
    cols = all_of(lipids_present),
    names_to = "lipid",
    values_to = "value"
  ) |>
  mutate(
    lipid = factor(lipid, levels = lipids_present),
    lipid_label = recode(as.character(lipid), !!!lipid_labels)
  )

sum_df <- df_long |>
  group_by(condition, infection, lipid, lipid_label) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = sd / sqrt(n),
    .groups = "drop"
  )

write.csv(
  sum_df,
  file = file.path(output_dir, "sms12_lipid_summary_statistics.csv"),
  row.names = FALSE
)

## ---- Plot helper ----
plot_one_lipid <- function(sum_df, lipid_name, y_label = "pmol / sample") {
  ggplot(
    sum_df |> filter(lipid == lipid_name),
    aes(x = condition, y = mean, fill = infection)
  ) +
    geom_col(
      position = position_dodge(width = 0.7),
      width = 0.6,
      colour = "black"
    ) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      position = position_dodge(width = 0.7),
      width = 0.2
    ) +
    scale_fill_manual(
      values = c("no" = "#009E73", "yes" = "#CC79A7"),
      labels = infection_labels,
      name = "Infection status"
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = y_label
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.key.size = unit(0.4, "cm")
    )
}

plots_by_lipid <- lapply(lipids_present, function(lipid_name) {
  plot_one_lipid(
    sum_df,
    lipid_name,
    y_label = unname(lipid_labels[[lipid_name]])
  )
})
names(plots_by_lipid) <- lipids_present

for (lipid_name in lipids_present) {
  ggsave(
    filename = file.path(output_dir, paste0(lipid_name, "_barplot.pdf")),
    plot = plots_by_lipid[[lipid_name]],
    width = 5,
    height = 4,
    units = "in"
  )

  ggsave(
    filename = file.path(output_dir, paste0(lipid_name, "_barplot.png")),
    plot = plots_by_lipid[[lipid_name]],
    width = 5,
    height = 4,
    units = "in",
    dpi = 300
  )
}

## ---- Combined faceted plot ----
p_combined <- ggplot(
  sum_df,
  aes(x = condition, y = mean, fill = infection)
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    colour = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  facet_wrap(~ lipid_label, scales = "free_y", nrow = 2) +
  scale_fill_manual(
    values = c("no" = "#009E73", "yes" = "#CC79A7"),
    labels = infection_labels,
    name = "Infection status"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05), add = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "pmol / sample"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "sms12_selected_lipid_barplots.pdf"),
  plot = p_combined,
  width = 11,
  height = 7,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "sms12_selected_lipid_barplots.png"),
  plot = p_combined,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300
)
