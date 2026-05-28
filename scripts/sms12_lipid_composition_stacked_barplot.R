library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ==============================
# SMS1/2 lipid composition stacked bar plot
# - Reads one raw Excel sheet containing infection, condition, and lipid columns
# - Averages selected lipids by condition and infection
# - Plots each lipid as a proportion of selected sphingolipid classes
# - Creates 100% stacked bar plots faceted by infection status
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
  "dh_sm_total"
)

lipid_labels <- c(
  "cer_total" = "Cer",
  "sm_total" = "SM",
  "hex_cer_total" = "HexCer",
  "lac_cer_total" = "LacCer",
  "dh_sm_total" = "dhSM"
)

lipid_palette <- c(
  "cer_total" = "#4B0092",
  "sm_total" = "#E66100",
  "hex_cer_total" = "#E1BE6A",
  "lac_cer_total" = "#40B0A6",
  "dh_sm_total" = "#D35FB7"
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

missing_lipid_labels <- setdiff(lipids_present, names(lipid_labels))
if (length(missing_lipid_labels) > 0) {
  warning(
    "These lipids do not have labels in lipid_labels; using column names: ",
    paste(missing_lipid_labels, collapse = ", ")
  )
  lipid_labels <- c(lipid_labels, setNames(missing_lipid_labels, missing_lipid_labels))
}

missing_lipid_colours <- setdiff(lipids_present, names(lipid_palette))
if (length(missing_lipid_colours) > 0) {
  warning(
    "These lipids do not have colours in lipid_palette; using fallback colours: ",
    paste(missing_lipid_colours, collapse = ", ")
  )
  lipid_palette <- c(
    lipid_palette,
    setNames(hue_pal()(length(missing_lipid_colours)), missing_lipid_colours)
  )
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
  infection_labels <- c(infection_labels, setNames(extra_infections, extra_infections))
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

## ---- Average selected lipids by condition and infection ----
df_mean_inf <- df_sms12_raw |>
  group_by(condition, infection) |>
  summarise(
    across(all_of(lipids_present), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(across(all_of(lipids_present), ~ ifelse(is.nan(.x), NA_real_, .x)))

write.csv(
  df_mean_inf,
  file = file.path(output_dir, "sms12_lipid_composition_averaged_values.csv"),
  row.names = FALSE
)

## ---- Reshape and calculate proportions ----
df_long_mean_inf <- df_mean_inf |>
  pivot_longer(
    cols = all_of(lipids_present),
    names_to = "lipid",
    values_to = "mean_value"
  ) |>
  mutate(
    lipid = factor(lipid, levels = lipids_present),
    lipid_label = recode(as.character(lipid), !!!lipid_labels),
    infection_label = recode(as.character(infection), !!!infection_labels),
    infection_label = factor(
      infection_label,
      levels = unname(infection_labels[infection_order])
    )
  ) |>
  group_by(condition, infection) |>
  mutate(
    total_selected_lipids = sum(mean_value, na.rm = TRUE),
    proportion = ifelse(total_selected_lipids > 0, mean_value / total_selected_lipids, NA_real_)
  ) |>
  ungroup()

write.csv(
  df_long_mean_inf,
  file = file.path(output_dir, "sms12_lipid_composition_long_values.csv"),
  row.names = FALSE
)

## ---- Stacked proportional bar plot ----
p_stack_inf <- ggplot(
  df_long_mean_inf,
  aes(x = condition, y = proportion, fill = lipid)
) +
  geom_col(position = "stack", colour = "black") +
  facet_wrap(~ infection_label, nrow = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = lipid_palette[lipids_present],
    labels = lipid_labels[lipids_present],
    name = "Sphingolipid class"
  ) +
  labs(
    x = NULL,
    y = "Proportion of sphingolipids"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  filename = file.path(output_dir, "sms12_lipid_composition_stacked_barplot.pdf"),
  plot = p_stack_inf,
  width = 9,
  height = 4.5,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "sms12_lipid_composition_stacked_barplot.png"),
  plot = p_stack_inf,
  width = 9,
  height = 4.5,
  units = "in",
  dpi = 300
)
