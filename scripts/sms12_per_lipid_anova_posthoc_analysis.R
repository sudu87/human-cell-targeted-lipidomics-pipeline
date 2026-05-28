library(readxl)
library(janitor)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(emmeans)
library(writexl)

# ==============================
# Per-lipid ANOVA and posthoc analysis
# - Reads one Excel sheet containing metadata and lipid abundance columns
# - Applies log10(x + 1) transformation to lipid columns
# - Runs per-lipid infection * condition models
# - Uses emmeans for posthoc contrasts and BH FDR correction
# ==============================

## ---- Configure input and outputs ----
# Edit these values before running the script.
input_file <- "path/to/input_lipidomics_data.xlsx"
sheet_name <- "sheet_name"

out_dir <- "analysis_outputs/statistics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Read and prepare data ----
sms_test <- read_excel(input_file, sheet = sheet_name) |>
  clean_names() |>
  mutate(
    infection = factor(infection, levels = c("no", "yes")),
    condition = factor(condition)
  )

meta_cols <- c("infection", "condition")
lipid_cols <- setdiff(names(sms_test), meta_cols)

sms_test <- sms_test |>
  mutate(across(all_of(lipid_cols), ~ suppressWarnings(as.numeric(.x)))) |>
  mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))

## ---- Per-lipid ANOVA table ----
per_lipid_anova <- function(lip) {
  fit <- lm(as.formula(paste(lip, "~ infection * condition")), data = sms_test)

  anova(fit) |>
    as.data.frame() |>
    rownames_to_column("term") |>
    as_tibble() |>
    mutate(lipid = lip, .before = 1)
}

anova_sms_test <- map_dfr(lipid_cols, per_lipid_anova)

write_xlsx(
  anova_sms_test,
  path = file.path(out_dir, "per_lipid_anova.xlsx")
)

## ---- Per-lipid posthoc contrasts ----
per_lipid_posthoc <- function(lip) {
  fit <- lm(as.formula(paste(lip, "~ infection * condition")), data = sms_test)

  condition_within_infection <- emmeans(fit, ~ condition | infection) |>
    contrast("pairwise", adjust = "none") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    mutate(lipid = lip, family = "Condition|within Infection", .before = 1)

  infection_within_condition <- emmeans(fit, ~ infection | condition) |>
    contrast("pairwise", adjust = "none") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    mutate(lipid = lip, family = "Infection|within Condition", .before = 1)

  out <- bind_rows(condition_within_infection, infection_within_condition)

  if (!"p.value" %in% names(out)) {
    if ("p_value" %in% names(out)) {
      out <- dplyr::rename(out, p.value = p_value)
    } else if (all(c("t.ratio", "df") %in% names(out))) {
      out$p.value <- 2 * pt(abs(out$`t.ratio`), df = out$df, lower.tail = FALSE)
    } else if ("z.ratio" %in% names(out)) {
      out$p.value <- 2 * pnorm(abs(out$`z.ratio`), lower.tail = FALSE)
    } else {
      stop("No p-values found; columns are: ", paste(names(out), collapse = ", "))
    }
  }

  out
}

res <- map_dfr(lipid_cols, per_lipid_posthoc)

## ---- BH correction across lipids for the same contrast within each stratum ----
cond_tbl <- res |>
  filter(family == "Condition|within Infection") |>
  drop_na(p.value) |>
  group_by(infection, contrast) |>
  mutate(FDR = p.adjust(p.value, "BH")) |>
  ungroup()

inf_tbl <- res |>
  filter(family == "Infection|within Condition") |>
  drop_na(p.value) |>
  group_by(condition, contrast) |>
  mutate(FDR = p.adjust(p.value, "BH")) |>
  ungroup()

posthoc_sms_test <- bind_rows(cond_tbl, inf_tbl)

write_xlsx(
  posthoc_sms_test,
  path = file.path(out_dir, "per_lipid_posthoc.xlsx")
)
