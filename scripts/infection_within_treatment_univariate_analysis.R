library(readxl)
library(janitor)
library(dplyr)
library(stringr)
library(purrr)
library(emmeans)
library(broom)
library(writexl)

## ---- Set your file paths here ----
file_path <- "path/to/your/input_file.xlsx"
output_file <- "path/to/your/output_file.xlsx"

## ---- Read the Excel file ----
df <- read_excel(file_path) %>%
  clean_names() %>%
  mutate(
    infection = str_to_title(as.character(infection)),
    infection = factor(infection, levels = c("Uninfected", "Infected"))
  )

## ---- Define lipid columns ----
lipid_cols <- setdiff(names(df), "infection")

## ---- Ensure numeric + log10 transform (BEFORE ANOVA) ----
df <- df %>%
  mutate(across(all_of(lipid_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(all_of(lipid_cols), ~ log10(. + 1)))

## ---- One-way ANOVA per lipid ----
get_infection_effect <- function(lip, data) {
  form <- reformulate("infection", response = lip)
  fit <- aov(form, data = data)

  em <- emmeans(fit, ~ infection)
  pw <- pairs(em, adjust = "fdr") %>% tidy()

  pw %>% mutate(lipid = lip, .before = 1)
}

## ---- Run across all lipids ----
res_infection <- map_dfr(lipid_cols, ~ get_infection_effect(.x, df))

head(res_infection)

## ---- Write results ----
write_xlsx(res_infection, output_file)
