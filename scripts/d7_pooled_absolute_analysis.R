# d7-dhSph labelled abundance: class pools and individual LCBs.
# Run from the repository root, or set ISOTOPE_ANALYSIS_ROOT.

required <- c("readxl", "dplyr", "tidyr", "emmeans", "ggplot2", "writexl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

library(dplyr)
library(tidyr)
library(ggplot2)

root <- normalizePath(Sys.getenv("ISOTOPE_ANALYSIS_ROOT", unset = getwd()), mustWork = TRUE)
absolute_file <- Sys.getenv(
  "D7_ABSOLUTE_LONG_FILE",
  unset = file.path(
    root, "outputs", "isotope_tracer_flux_analysis", "flux_absolute_long.csv"
  )
)
workbook_dir <- Sys.getenv("D7_POOLED_WORKBOOK_DIR", unset = file.path(root, "data", "validation", "d7_dhsph"))
output_dir <- Sys.getenv(
  "D7_POOLED_OUTPUT_DIR",
  unset = file.path(
    root, "outputs", "isotope_tracer_absolute_labelled_abundance", "d7_dhsph"
  )
)
pseudovalue <- as.numeric(Sys.getenv("D7_LOG10_PSEUDOVALUE", unset = "0.01"))
if (!file.exists(absolute_file)) stop("Missing absolute input CSV: ", absolute_file)
if (!is.finite(pseudovalue) || pseudovalue <= 0) stop("Pseudovalue must be positive")

classes <- c("dhCer", "Cer", "dhSM", "SM", "HexCer", "LacCer")
stems <- setNames(c("dhcer", "cer", "dhsm", "sm", "hexcer", "laccer"), classes)
lcb_species <- c("d7-dhSph", "d7-Sph", "d7-S1P")
component_species <- c("d7-16:0", "d7-24:0", "d7-24:1")
times <- c(0.5, 1, 2, 4, 8)

absolute <- read.csv(absolute_file, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("experiment", "unit", "berlin_sample_id", "infection", "is_tracer_treated",
            "time_h", "replicate", "lipid_class_short", "species", "value", "is_class_total")
if (length(setdiff(needed, names(absolute)))) stop("Absolute CSV lacks required columns")

d7 <- absolute |>
  filter(experiment == "d7-dhSph", unit == "pmol_per_sample", is_tracer_treated,
         !is_class_total, infection %in% c("Uninfected", "Infected"), time_h %in% times) |>
  mutate(replicate = as.character(replicate))

components <- d7 |>
  filter(lipid_class_short %in% classes, species %in% component_species) |>
  select(berlin_sample_id, infection, time_h, replicate, lipid_class_short, species, value) |>
  pivot_wider(names_from = species, values_from = value)
if (nrow(components) != 6 * 30 || anyNA(components[component_species])) {
  stop("Expected all three d7 components for every class and sample")
}

pooled <- components |>
  transmute(berlin_sample_id, infection, time_h, replicate, lipid_class_short,
            feature_type = "Pooled labelled class", feature_name = paste0("d7-", lipid_class_short),
            component_d7_16_0 = `d7-16:0`, component_d7_24_0 = `d7-24:0`,
            component_d7_24_1 = `d7-24:1`,
            value = component_d7_16_0 + component_d7_24_0 + component_d7_24_1)

lcb <- d7 |>
  filter(lipid_class_short == "LCB", species %in% lcb_species) |>
  transmute(berlin_sample_id, infection, time_h, replicate, lipid_class_short,
            feature_type = "Individual labelled LCB", feature_name = species,
            component_d7_16_0 = NA_real_, component_d7_24_0 = NA_real_,
            component_d7_24_1 = NA_real_, value)

readouts <- bind_rows(lcb, pooled) |>
  mutate(experiment = "d7-dhSph", unit = "pmol_per_sample", pseudovalue = pseudovalue,
         value_log10 = log10(value + pseudovalue)) |>
  select(experiment, berlin_sample_id, infection, time_h, replicate, feature_type,
         lipid_class_short, feature_name, unit, starts_with("component_"), value,
         pseudovalue, value_log10) |>
  arrange(feature_name, time_h, infection, replicate)

design <- readouts |> count(feature_name, time_h, infection, name = "n")
if (nrow(readouts) != 270 || nrow(design) != 90 || any(design$n != 3) ||
    any(!is.finite(readouts$value)) || any(readouts$value < 0)) {
  stop("Expected nine endpoints, five times, two infections, n=3; all values finite and nonnegative")
}

check_workbook <- function(lipid_class) {
  file <- file.path(workbook_dir, paste0("d7_dhsph_", stems[[lipid_class]],
                                          "_d7_sum_16_0_24_0_24_1.xlsx"))
  if (!file.exists(file)) stop("Missing workbook: ", file)
  x <- readxl::read_excel(file, sheet = "Raw d7 species", range = "A1:G31")
  names(x) <- c("time_h", "infection_label", "replicate", "wb_16", "wb_24_0", "wb_24_1", "wb_total")
  x <- x |>
    mutate(infection = recode(infection_label, `-Sne` = "Uninfected", `+Sne` = "Infected"),
           replicate = as.character(replicate))
  source <- components |>
    filter(lipid_class_short == lipid_class) |>
    select(time_h, infection, replicate, `d7-16:0`, `d7-24:0`, `d7-24:1`)
  joined <- x |> left_join(source, by = c("time_h", "infection", "replicate"))
  if (nrow(joined) != 30 || anyNA(joined[, c("wb_16", "wb_24_0", "wb_24_1", "wb_total",
                                          "d7-16:0", "d7-24:0", "d7-24:1")])) {
    stop("Incomplete or mismatched workbook: ", basename(file))
  }
  sum_delta <- max(abs(joined$wb_total - joined$wb_16 - joined$wb_24_0 - joined$wb_24_1))
  source_delta <- max(abs(joined$wb_16 - joined$`d7-16:0`),
                      abs(joined$wb_24_0 - joined$`d7-24:0`),
                      abs(joined$wb_24_1 - joined$`d7-24:1`),
                      abs(joined$wb_total - joined$`d7-16:0` - joined$`d7-24:0` - joined$`d7-24:1`))
  if (sum_delta > 1e-9 || source_delta > 1e-9) stop("Source mismatch: ", basename(file))
  tibble(feature_name = paste0("d7-", lipid_class), workbook = basename(file),
         n = nrow(joined), zeros = sum(joined$wb_total == 0),
         max_sum_delta = sum_delta, max_source_delta = source_delta)
}

expected_validation_files <- c(
  file.path(workbook_dir, paste0("d7_dhsph_", unname(stems), "_d7_sum_16_0_24_0_24_1.xlsx")),
  file.path(workbook_dir, paste0("d7_dhsph_", c("d7_dhsph", "d7_sph", "d7_s1p"), ".xlsx"))
)
available_validation_files <- file.exists(expected_validation_files)
if (any(available_validation_files) && !all(available_validation_files)) {
  stop("d7 validation directory is incomplete; provide all nine workbooks or none.")
}
run_workbook_qc <- all(available_validation_files)
qc <- if (run_workbook_qc) bind_rows(lapply(classes, check_workbook)) else tibble(
  feature_name = paste0("d7-", classes), workbook = NA_character_, n = NA_integer_,
  zeros = NA_integer_, max_sum_delta = NA_real_, max_source_delta = NA_real_,
  note = "Optional validation workbooks not supplied; primary source-table checks still ran."
)

prism_files <- c(
  setNames(paste0("d7_dhsph_", c("d7_dhsph", "d7_sph", "d7_s1p"), ".xlsx"), lcb_species),
  setNames(paste0("d7_dhsph_", unname(stems), "_d7_sum_16_0_24_0_24_1.xlsx"),
           paste0("d7-", classes))
)
check_prism <- function(name) {
  file <- file.path(workbook_dir, prism_files[[name]])
  x <- readxl::read_excel(file, sheet = 1, range = "A1:G6", col_names = FALSE,
                          .name_repair = "minimal")
  if (nrow(x) != 6 || ncol(x) != 7 || !identical(as.character(unlist(x[1, c(2, 5)])),
                                                  c("-Sne", "+Sne"))) {
    stop("Unexpected Prism layout: ", basename(file))
  }
  values <- bind_rows(lapply(seq_along(times), function(i) {
    bind_rows(lapply(2:7, function(j) {
      tibble(time_h = times[i], infection = if (j <= 4) "Uninfected" else "Infected",
             replicate = as.character(if (j <= 4) j - 1 else j - 4),
             prism_value = as.numeric(x[[j]][i + 1]))
    }))
  }))
  source <- readouts |>
    filter(feature_name == name) |>
    select(time_h, infection, replicate, value)
  joined <- values |> left_join(source, by = c("time_h", "infection", "replicate"))
  if (nrow(joined) != 30 || anyNA(joined) || max(abs(joined$prism_value - joined$value)) > 1e-9) {
    stop("Prism/source mismatch: ", basename(file))
  }
  tibble(feature_name = name, workbook = basename(file), n = nrow(joined),
         max_prism_source_delta = max(abs(joined$prism_value - joined$value)))
}
prism_qc <- if (run_workbook_qc) bind_rows(lapply(names(prism_files), check_prism)) else tibble(
  feature_name = names(prism_files), workbook = NA_character_, n = NA_integer_,
  max_prism_source_delta = NA_real_,
  note = "Optional Prism validation workbooks not supplied."
)

fit_model <- function(data, offset = pseudovalue) {
  data <- data |>
    mutate(infection = factor(infection, levels = c("Uninfected", "Infected")),
           time_f = factor(time_h, levels = times), response = log10(value + offset))
  lm(response ~ infection * time_f, data = data)
}

pairwise_for_feature <- function(data) {
  fit <- fit_model(data)
  em <- emmeans::emmeans(fit, ~ infection | time_f) |>
    emmeans::contrast(method = list("Infected - Uninfected" = c(-1, 1)), adjust = "none") |>
    summary(infer = c(TRUE, TRUE)) |>
    as_tibble()
  observed <- data |>
    group_by(time_h) |>
    summarise(n_uninfected = sum(infection == "Uninfected"), n_infected = sum(infection == "Infected"),
              mean_uninfected_pmol = mean(value[infection == "Uninfected"]),
              mean_infected_pmol = mean(value[infection == "Infected"]),
              all_observed_zero = all(value == 0), .groups = "drop")
  em |>
    transmute(time_h = as.numeric(as.character(time_f)), contrast = as.character(contrast),
              estimate_log10_infected_minus_uninfected = estimate, se_log10 = SE, df,
              t_ratio = t.ratio, lower_cl_log10 = lower.CL, upper_cl_log10 = upper.CL,
              ratio_infected_over_uninfected = 10^estimate,
              ratio_lower_95 = 10^lower.CL, ratio_upper_95 = 10^upper.CL,
              p_value = p.value) |>
    left_join(observed, by = "time_h") |>
    mutate(test_status = if_else(all_observed_zero, "all_zero_observed_model_based", "estimated"),
           model = "log10(value + pseudovalue) ~ infection * categorical_time") |>
    select(-all_observed_zero)
}

pairwise_raw <- readouts |>
  group_by(experiment, feature_type, lipid_class_short, feature_name) |>
  group_modify(~ pairwise_for_feature(.x)) |>
  ungroup() |>
  arrange(feature_name, time_h)
if (nrow(pairwise_raw) != 45 || any(!is.finite(pairwise_raw$p_value))) stop("Expected 45 finite contrasts")

pairwise_fdr <- pairwise_raw |>
  mutate(p_fdr = p.adjust(p_value, method = "BH"),
         fdr_method = "BH across 45 d7 endpoint/time contrasts",
         significant_fdr_0_05 = p_fdr < 0.05)

curve_for_feature <- function(data) {
  fit <- fit_model(data)
  a <- anova(fit)
  tibble(n_samples = nrow(data), model = "log10(value + pseudovalue) ~ infection * categorical_time",
         p_infection = a["infection", "Pr(>F)"], p_time = a["time_f", "Pr(>F)"],
         p_infection_time = a["infection:time_f", "Pr(>F)"],
         shapiro_residual_p = shapiro.test(residuals(fit))$p.value,
         max_abs_standardized_residual = max(abs(rstandard(fit)), na.rm = TRUE),
         max_cooks_distance = max(cooks.distance(fit), na.rm = TRUE))
}

curves <- readouts |>
  group_by(experiment, feature_type, lipid_class_short, feature_name) |>
  group_modify(~ curve_for_feature(.x)) |>
  ungroup() |>
  mutate(p_fdr_infection = p.adjust(p_infection, method = "BH"),
         p_fdr_time = p.adjust(p_time, method = "BH"),
         curve_fdr = p.adjust(p_infection_time, method = "BH"),
         significant_infection_fdr_0_05 = p_fdr_infection < 0.05,
         significant_time_fdr_0_05 = p_fdr_time < 0.05,
         significant_curve_fdr_0_05 = curve_fdr < 0.05,
         fdr_method = "Separate BH families of nine endpoints for infection, time, and interaction") |>
  arrange(feature_name)

curve_labels <- curves |>
  transmute(feature_name, curve_fdr,
            curve_fdr_label = if_else(curve_fdr < 0.001, "FDR < 0.001",
                                      sprintf("FDR = %.3f", curve_fdr)))

sensitivity <- readouts |>
  filter(feature_name == "d7-LacCer") |>
  group_by(feature_name) |>
  group_modify(~ bind_rows(lapply(c(0.001, 0.01, 0.1), function(offset) {
    fit <- fit_model(.x, offset)
    em <- emmeans::emmeans(fit, ~ infection | time_f) |>
      emmeans::contrast(method = list("Infected - Uninfected" = c(-1, 1)), adjust = "none") |>
      summary(infer = c(TRUE, TRUE)) |>
      as_tibble()
    tibble(pseudovalue = offset,
           p_infection_time = anova(fit)["infection:time_f", "Pr(>F)"],
           p_2h = em$p.value[as.character(em$time_f) == "2"])
  }))) |>
  ungroup()

summary_values <- readouts |>
  group_by(feature_name, time_h, infection) |>
  summarise(mean = mean(value), sd = sd(value), n = n(), .groups = "drop") |>
  left_join(curve_labels, by = "feature_name") |>
  mutate(infection = factor(infection, levels = c("Uninfected", "Infected")),
         time_f = factor(time_h, levels = times),
         panel = paste(feature_name, curve_fdr_label, sep = "\n"))

plot_base <- function(data) {
  ggplot(data, aes(time_f, mean, group = infection, colour = infection, shape = infection)) +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd), width = 0.12, linewidth = 0.7) +
    geom_point(size = 3) +
    scale_colour_manual(values = c(Uninfected = "#0000FF", Infected = "#FF0000"),
                        labels = c(Uninfected = "-Sne", Infected = "+Sne")) +
    scale_shape_manual(values = c(Uninfected = 16, Infected = 15),
                       labels = c(Uninfected = "-Sne", Infected = "+Sne")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Time (h)", y = "Labelled abundance (pmol/sample)", colour = NULL, shape = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(output_dir, "figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(readouts, file.path(output_dir, "d7_pooled_absolute_readouts.csv"), row.names = FALSE)
write.csv(qc, file.path(output_dir, "d7_workbook_qc.csv"), row.names = FALSE)
write.csv(prism_qc, file.path(output_dir, "d7_prism_crosscheck_qc.csv"), row.names = FALSE)
write.csv(pairwise_raw, file.path(output_dir, "timepoint_pairwise_raw.csv"), row.names = FALSE)
write.csv(pairwise_fdr, file.path(output_dir, "timepoint_pairwise_experiment_fdr.csv"), row.names = FALSE)
writexl::write_xlsx(pairwise_fdr, file.path(output_dir, "timepoint_pairwise_experiment_fdr.xlsx"))
write.csv(curves, file.path(output_dir, "infection_time_curve_models.csv"), row.names = FALSE)
writexl::write_xlsx(curves, file.path(output_dir, "infection_time_curve_models.xlsx"))
write.csv(curve_labels, file.path(output_dir, "curve_fdr_labels.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(output_dir, "laccer_pseudovalue_sensitivity.csv"), row.names = FALSE)
write.csv(summary_values, file.path(output_dir, "plot_summary_mean_sd.csv"), row.names = FALSE)

for (name in unique(summary_values$feature_name)) {
  p <- plot_base(filter(summary_values, feature_name == name)) +
    labs(title = paste0(name, " labelled abundance"))
  ggsave(file.path(plot_dir, paste0(gsub("[^A-Za-z0-9]+", "_", name), ".pdf")),
         p, width = 6.4, height = 4.4)
}
combined <- plot_base(summary_values) + facet_wrap(~ panel, scales = "free_y", ncol = 3)
ggsave(file.path(plot_dir, "d7_pooled_timecourses_combined.pdf"), combined,
       width = 13, height = 12)

pdf(file.path(output_dir, "model_residual_diagnostics.pdf"), width = 9, height = 4)
for (name in unique(readouts$feature_name)) {
  data <- filter(readouts, feature_name == name)
  fit <- fit_model(data)
  par(mfrow = c(1, 2))
  plot(fitted(fit), residuals(fit), xlab = "Fitted", ylab = "Residual", main = name)
  abline(h = 0, lty = 2)
  qqnorm(residuals(fit), main = paste(name, "Q-Q"))
  qqline(residuals(fit))
}
dev.off()

message("Output directory: ", output_dir)
message("Readouts: ", nrow(readouts), "; pairwise tests: ", nrow(pairwise_fdr))
message("Timepoint FDR < 0.05: ", sum(pairwise_fdr$significant_fdr_0_05))
message("Curve FDR < 0.05: ", sum(curves$significant_curve_fdr_0_05))
