library(readxl)
library(janitor)
library(dplyr)
library(purrr)
library(emmeans)
library(broom)
library(writexl)
library(rcompanion)

# ==============================
# Simple total lipid MANOVA script
# - Reads one Excel file
# - Runs MANOVA using Pillai's trace
# - Runs per-lipid ANOVA
# - Runs posthoc contrasts with BH/FDR correction
# ==============================

## ---- Set input and output locations here ----
input_file <- "path/to/input_file.xlsx"
out_dir <- "path/to/output_directory"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_infection <- file.path(
  out_dir,
  "posthoc_infection_contrasts_by_inhibitor.xlsx"
)

out_inhibitor <- file.path(
  out_dir,
  "posthoc_inhibitor_contrasts_by_infection.xlsx"
)

## ---- Read the Excel file ----
df <- read_excel(input_file) %>%
  clean_names() %>%
  mutate(
    infection = factor(infection),
    inhibitor = factor(inhibitor)
  )

## ---- Define lipid columns ----
meta_cols  <- c("infection", "inhibitor")
lipid_cols <- setdiff(names(df), meta_cols)

## ---- Log10 transform lipid values ----
df <- df %>%
  mutate(across(all_of(lipid_cols), ~ log10(as.numeric(.x) + 1)))

## ---- MANOVA ----
resp_mat <- as.matrix(df[, lipid_cols])

manova_fit <- manova(resp_mat ~ infection * inhibitor, data = df)

summary(manova_fit, test = "Pillai")

summary.aov(manova_fit)

## ---- Simple residual normality plot ----
plotNormalDensity(as.vector(residuals(manova_fit)))

## ---- Per-lipid ANOVA: infection contrasts within each inhibitor ----
get_infection_within_treatment <- function(lip, data) {
  form <- as.formula(paste(lip, "~ infection * inhibitor"))
  fit  <- aov(form, data = data)

  em <- emmeans(fit, ~ infection | inhibitor)
  pw <- pairs(em, adjust = "BH") %>% tidy()

  pw %>% mutate(lipid = lip, .before = 1)
}

res_infection <- map_dfr(
  lipid_cols,
  ~ get_infection_within_treatment(.x, df)
)

write_xlsx(res_infection, out_infection)

## ---- Per-lipid ANOVA: inhibitor contrasts within each infection group ----
get_treatment_within_infection <- function(lip, data) {
  form <- as.formula(paste(lip, "~ infection * inhibitor"))
  fit  <- aov(form, data = data)

  em <- emmeans(fit, ~ inhibitor | infection)
  pw <- pairs(em, adjust = "BH") %>% tidy()

  pw %>% mutate(lipid = lip, .before = 1)
}

res_inhibitor <- map_dfr(
  lipid_cols,
  ~ get_treatment_within_infection(.x, df)
)

write_xlsx(res_inhibitor, out_inhibitor)
