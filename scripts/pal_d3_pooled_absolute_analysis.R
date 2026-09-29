# Pal-d3 16:0 d3+d6 labelled abundance: infected vs uninfected by timepoint.
# Run from the repository root, or set ISOTOPE_ANALYSIS_ROOT.

required_packages <- c("readxl", "dplyr", "tidyr", "emmeans", "writexl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(tidyr)

analysis_root <- normalizePath(
  Sys.getenv("ISOTOPE_ANALYSIS_ROOT", unset = getwd()),
  mustWork = TRUE
)
absolute_file <- Sys.getenv(
  "PAL_D3_ABSOLUTE_LONG_FILE",
  unset = file.path(
    analysis_root, "outputs", "isotope_tracer_flux_analysis",
    "flux_absolute_long.csv"
  )
)
combined_workbook_dir <- Sys.getenv(
  "PAL_D3_COMBINED_WORKBOOK_DIR",
  unset = file.path(analysis_root, "data", "validation", "pal_d3")
)
old_pairwise_file <- Sys.getenv(
  "PAL_D3_OLD_PAIRWISE_FILE",
  unset = file.path(analysis_root, "outputs", "isotope_tracer_flux_analysis",
                    "timepoint_pairwise_emmeans_experiment_fdr_labelled_incorporation.csv")
)
output_dir <- Sys.getenv(
  "PAL_D3_D3_D6_OUTPUT_DIR",
  unset = file.path(
    analysis_root, "outputs", "isotope_tracer_absolute_labelled_abundance",
    "pal_d3"
  )
)
pseudovalue <- as.numeric(Sys.getenv("PAL_D3_LOG10_PSEUDOVALUE", unset = "0.01"))

if (!file.exists(absolute_file)) stop("Missing absolute input CSV: ", absolute_file)
if (!is.finite(pseudovalue) || pseudovalue <= 0) {
  stop("PAL_D3_LOG10_PSEUDOVALUE must be positive.")
}

class_names <- c("dhCer", "Cer", "dhSM", "SM", "HexCer", "LacCer")
class_file_stems <- c("dhcer", "cer", "dhsm", "sm", "hexcer", "laccer")
names(class_file_stems) <- class_names
timepoints <- c(0.5, 1, 2, 4, 8)

absolute <- read.csv(absolute_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c(
  "experiment", "unit", "berlin_sample_id", "infection", "is_tracer_treated",
  "time_h", "replicate", "lipid_class_short", "species", "value", "is_class_total"
)
missing_cols <- setdiff(required_cols, names(absolute))
if (length(missing_cols)) {
  stop("Missing columns in absolute input: ", paste(missing_cols, collapse = ", "))
}

pal <- absolute |>
  filter(
    experiment == "Pal-d3",
    unit == "pmol_per_sample",
    is_tracer_treated,
    !is_class_total,
    infection %in% c("Uninfected", "Infected"),
    time_h %in% timepoints
  ) |>
  mutate(replicate = as.character(replicate))

paired <- pal |>
  filter(
    lipid_class_short %in% class_names,
    species %in% c("d3-16:0", "d6-16:0")
  ) |>
  select(
    berlin_sample_id, infection, time_h, replicate, lipid_class_short,
    species, value
  ) |>
  pivot_wider(names_from = species, values_from = value)

if (anyNA(paired$`d3-16:0`) || anyNA(paired$`d6-16:0`)) {
  stop("A paired d3/d6 measurement is missing in the absolute source table.")
}

paired <- paired |>
  transmute(
    berlin_sample_id, infection, time_h, replicate, lipid_class_short,
    feature_type = "Lipid species",
    feature_name = paste(lipid_class_short, "d3+d6-16:0"),
    component_d3 = `d3-16:0`,
    component_d6 = `d6-16:0`,
    value = component_d3 + component_d6
  )

lcb <- pal |>
  filter(lipid_class_short == "LCB", species %in% c("d3-dhSph", "d3-Sph")) |>
  transmute(
    berlin_sample_id, infection, time_h, replicate, lipid_class_short,
    feature_type = "LCB molecule",
    feature_name = species,
    component_d3 = value,
    component_d6 = NA_real_,
    value
  )

readouts <- bind_rows(lcb, paired) |>
  mutate(
    experiment = "Pal-d3",
    unit = "pmol_per_sample",
    pseudovalue = pseudovalue,
    value_log10 = log10(value + pseudovalue)
  ) |>
  select(
    experiment, berlin_sample_id, infection, time_h, replicate,
    feature_type, lipid_class_short, feature_name, unit,
    component_d3, component_d6, value, pseudovalue, value_log10
  ) |>
  arrange(feature_name, time_h, infection, replicate)

if (any(!is.finite(readouts$value)) || any(readouts$value < 0)) {
  stop("Primary readouts must be finite and non-negative.")
}

design <- readouts |>
  count(feature_name, time_h, infection, name = "n")
if (
  n_distinct(readouts$feature_name) != 8 ||
  nrow(readouts) != 8 * 5 * 2 * 3 ||
  nrow(design) != 8 * 5 * 2 ||
  any(design$n != 3)
) {
  stop("The primary readouts do not match eight features, five times, and n=3 per group.")
}

expected_validation_files <- file.path(
  combined_workbook_dir,
  paste0("pal_d3_", unname(class_file_stems), "_d3_plus_d6_16_0.xlsx")
)
available_validation_files <- file.exists(expected_validation_files)
if (any(available_validation_files) && !all(available_validation_files)) {
  stop("Pal-d3 validation directory is incomplete; provide all six workbooks or none.")
}
run_workbook_qc <- all(available_validation_files)
qc <- if (run_workbook_qc) lapply(class_names, function(lipid_class) {
  file <- file.path(
    combined_workbook_dir,
    paste0("pal_d3_", class_file_stems[[lipid_class]], "_d3_plus_d6_16_0.xlsx")
  )
  if (!file.exists(file)) stop("Missing combined workbook: ", file)

  workbook <- readxl::read_excel(file, sheet = "Raw d3+d6", range = "A1:F31") |>
    rename(
      time_h = `Time (h)`,
      infection_label = Infection,
      replicate = Replicate,
      workbook_d3 = `d3 value`,
      workbook_d6 = `d6 value`,
      workbook_sum = `d3+d6 value`
    ) |>
    mutate(
      infection = recode(infection_label, `-Sne` = "Uninfected", `+Sne` = "Infected"),
      replicate = as.character(replicate)
    )

  if (nrow(workbook) != 30 || anyNA(workbook[, 1:6])) {
    stop("Incomplete raw sheet in ", basename(file))
  }

  source <- paired |>
    filter(lipid_class_short == lipid_class) |>
    select(
      time_h, infection, replicate,
      source_d3 = component_d3,
      source_d6 = component_d6,
      source_sum = value
    )
  checked <- workbook |>
    left_join(source, by = c("time_h", "infection", "replicate"))

  if (nrow(checked) != 30 || anyNA(checked$source_sum)) {
    stop("Workbook rows do not match source samples in ", basename(file))
  }

  sum_delta <- max(abs(checked$workbook_sum - checked$workbook_d3 - checked$workbook_d6))
  source_delta <- max(
    abs(checked$workbook_d3 - checked$source_d3),
    abs(checked$workbook_d6 - checked$source_d6),
    abs(checked$workbook_sum - checked$source_sum)
  )
  if (sum_delta > 1e-9 || source_delta > 1e-9) {
    stop("d3+d6 or source-table mismatch in ", basename(file))
  }

  tibble(
    lipid_class_short = lipid_class,
    workbook = basename(file),
    n_rows = nrow(checked),
    n_zero_sums = sum(checked$workbook_sum == 0),
    max_sum_delta = sum_delta,
    max_source_delta = source_delta
  )
}) |> bind_rows() else tibble(
  lipid_class_short = class_names,
  workbook = NA_character_, n_rows = NA_integer_, n_zero_sums = NA_integer_,
  max_sum_delta = NA_real_, max_source_delta = NA_real_,
  note = "Optional validation workbooks not supplied; primary source-table checks still ran."
)

fit_interaction_model <- function(data) {
  data <- data |>
    mutate(
      infection = factor(infection, levels = c("Uninfected", "Infected")),
      time_f = factor(time_h, levels = timepoints)
    )
  lm(value_log10 ~ infection * time_f, data = data)
}

fit_timepoint_contrasts <- function(data) {
  fit <- fit_interaction_model(data)
  contrasts <- emmeans::emmeans(fit, ~ infection | time_f) |>
    emmeans::contrast(
      method = list("Infected - Uninfected" = c(-1, 1)),
      adjust = "none"
    ) |>
    summary(infer = c(TRUE, TRUE)) |>
    as_tibble()

  observed <- data |>
    group_by(time_h) |>
    summarise(
      n_uninfected = sum(infection == "Uninfected"),
      n_infected = sum(infection == "Infected"),
      mean_uninfected_pmol = mean(value[infection == "Uninfected"]),
      mean_infected_pmol = mean(value[infection == "Infected"]),
      all_observed_zero = all(value == 0),
      .groups = "drop"
    )

  contrasts |>
    transmute(
      time_h = as.numeric(as.character(time_f)),
      contrast = as.character(contrast),
      estimate_log10_infected_minus_uninfected = estimate,
      se_log10 = SE,
      df,
      t_ratio = t.ratio,
      lower_cl_log10 = lower.CL,
      upper_cl_log10 = upper.CL,
      ratio_infected_over_uninfected = 10^estimate,
      ratio_lower_95 = 10^lower.CL,
      ratio_upper_95 = 10^upper.CL,
      p_value = p.value
    ) |>
    left_join(observed, by = "time_h") |>
    mutate(
      test_status = if_else(
        all_observed_zero,
        "all_zero_observed_model_based",
        "estimated"
      ),
      model = "log10(value + pseudovalue) ~ infection * categorical_time"
    ) |>
    select(-all_observed_zero)
}

pairwise_raw <- readouts |>
  group_by(experiment, feature_type, lipid_class_short, feature_name) |>
  group_modify(~ fit_timepoint_contrasts(.x)) |>
  ungroup() |>
  arrange(feature_name, time_h)

if (nrow(pairwise_raw) != 40 || any(!is.finite(pairwise_raw$p_value))) {
  stop("Expected 40 finite interaction-model timepoint contrasts.")
}

pairwise_fdr <- pairwise_raw |>
  mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    fdr_method = "BH across 40 Pal-d3 d3+d6 endpoint/time contrasts",
    significant_fdr_0_05 = p_fdr < 0.05
  )

fit_curve <- function(data) {
  fit <- fit_interaction_model(data)
  model_anova <- anova(fit)
  tibble(
    n_samples = nrow(data),
    model = "log10(value + pseudovalue) ~ infection * categorical_time",
    p_infection = model_anova["infection", "Pr(>F)"],
    p_time = model_anova["time_f", "Pr(>F)"],
    p_infection_time = model_anova["infection:time_f", "Pr(>F)"]
  )
}

curve_models <- readouts |>
  group_by(experiment, feature_type, lipid_class_short, feature_name) |>
  group_modify(~ fit_curve(.x)) |>
  ungroup() |>
  mutate(
    p_fdr_infection = p.adjust(p_infection, method = "BH"),
    p_fdr_time = p.adjust(p_time, method = "BH"),
    p_fdr_infection_time = p.adjust(p_infection_time, method = "BH"),
    significant_curve_fdr_0_05 = p_fdr_infection_time < 0.05,
    fdr_method = "BH across eight Pal-d3 pooled absolute endpoints"
  ) |>
  arrange(feature_name)

format_curve_fdr <- function(p) {
  ifelse(
    is.na(p),
    "curve FDR = NA",
    ifelse(p < 0.001, "curve FDR < 0.001", sprintf("curve FDR = %.3f", p))
  )
}

curve_labels <- curve_models |>
  transmute(
    experiment, feature_type, feature_name,
    curve_fdr = p_fdr_infection_time,
    curve_fdr_label = format_curve_fdr(p_fdr_infection_time)
  )

old_comparison <- tibble()
if (file.exists(old_pairwise_file)) {
  old <- read.csv(old_pairwise_file, stringsAsFactors = FALSE) |>
    filter(experiment == "Pal-d3") |>
    mutate(
      pooled_feature_name = if_else(
        feature_type == "LCB molecule",
        feature_name,
        sub(" d[36]-16:0$", " d3+d6-16:0", feature_name)
      ),
      old_component = case_when(
        feature_type == "LCB molecule" ~ "d3 LCB",
        grepl(" d3-16:0$", feature_name) ~ "d3 species",
        TRUE ~ "d6 species"
      )
    ) |>
    select(
      pooled_feature_name, time_h, old_feature_name = feature_name,
      old_component, old_p_value = p_value, old_p_fdr = p_fdr
    )

  old_comparison <- pairwise_fdr |>
    select(
      pooled_feature_name = feature_name, time_h,
      pooled_p_value = p_value, pooled_p_fdr = p_fdr,
      pooled_test_status = test_status
    ) |>
    left_join(old, by = c("pooled_feature_name", "time_h")) |>
    arrange(pooled_feature_name, time_h, old_component)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(readouts, file.path(output_dir, "pal_d3_d3_plus_d6_absolute_readouts.csv"), row.names = FALSE)
write.csv(qc, file.path(output_dir, "combined_workbook_qc.csv"), row.names = FALSE)
write.csv(pairwise_raw, file.path(output_dir, "timepoint_pairwise_emmeans_raw.csv"), row.names = FALSE)
write.csv(pairwise_fdr, file.path(output_dir, "timepoint_pairwise_emmeans_experiment_fdr.csv"), row.names = FALSE)
writexl::write_xlsx(pairwise_fdr, file.path(output_dir, "timepoint_pairwise_emmeans_experiment_fdr.xlsx"))
write.csv(curve_models, file.path(output_dir, "infection_time_curve_models.csv"), row.names = FALSE)
writexl::write_xlsx(curve_models, file.path(output_dir, "infection_time_curve_models.xlsx"))
write.csv(curve_labels, file.path(output_dir, "curve_fdr_labels.csv"), row.names = FALSE)
if (nrow(old_comparison)) {
  write.csv(old_comparison, file.path(output_dir, "old_vs_pooled_timepoint_comparison.csv"), row.names = FALSE)
}

message("Output directory: ", output_dir)
message("Readouts: ", nrow(readouts), "; planned tests: ", nrow(pairwise_fdr))
message("Estimable tests: ", sum(is.finite(pairwise_fdr$p_value)))
message("FDR < 0.05: ", sum(pairwise_fdr$significant_fdr_0_05))
message("Curve FDR < 0.05: ", sum(curve_models$significant_curve_fdr_0_05))
