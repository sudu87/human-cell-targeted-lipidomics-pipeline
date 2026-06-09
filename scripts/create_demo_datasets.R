library(writexl)

# ==============================
# Create small demo lipidomics datasets
# - Generates deterministic SMS1/2-style demo data
# - Writes one Excel workbook with raw total lipids, raw individual lipids,
#   and pairwise contrast tables
# - Intended for installation/demo checks, not biological interpretation
# ==============================

set.seed(42)

out_dir <- "demo_data"
out_file <- file.path(out_dir, "sms12_demo_lipidomics.xlsx")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

conditions <- c(
  "sms12_WT",
  "sms12_DKO",
  "sms1_16_no_AHT",
  "sms1_16_AHT",
  "sms2_no_AHT",
  "sms2_AHT",
  "sms2-M64R_no_AHT",
  "sms2-M64R_AHT"
)

infections <- c("no", "yes")
replicates <- 1:3

design <- expand.grid(
  condition = conditions,
  infection = infections,
  replicate = replicates,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

design <- design[
  order(
    match(design$condition, conditions),
    match(design$infection, infections),
    design$replicate
  ),
]

condition_scale <- c(
  "sms12_WT" = 1.00,
  "sms12_DKO" = 0.72,
  "sms1_16_no_AHT" = 0.86,
  "sms1_16_AHT" = 1.18,
  "sms2_no_AHT" = 0.90,
  "sms2_AHT" = 1.12,
  "sms2-M64R_no_AHT" = 0.80,
  "sms2-M64R_AHT" = 0.98
)

infection_scale <- c("no" = 1.00, "yes" = 1.25)

make_measurement <- function(base, lipid_shift = 0) {
  scale <- condition_scale[design$condition] * infection_scale[design$infection]
  replicate_noise <- c(-0.05, 0.01, 0.06)[match(design$replicate, replicates)]
  random_noise <- rnorm(nrow(design), mean = 0, sd = 0.025)
  value <- base * scale * (1 + lipid_shift + replicate_noise + random_noise)
  round(pmax(value, 0), 3)
}

total_specs <- c(
  "dh_sph" = 0.65,
  "sph" = 1.20,
  "s1p" = 0.42,
  "dh_cer_total" = 2.50,
  "cer_total" = 6.20,
  "dh_sm_total" = 3.15,
  "sm_total" = 18.50,
  "hex_cer_total" = 4.80,
  "lac_cer_total" = 2.10
)

sms12_total_demo <- design[c("infection", "condition")]
for (lipid in names(total_specs)) {
  shift <- (match(lipid, names(total_specs)) - 5) * 0.015
  sms12_total_demo[[lipid]] <- make_measurement(total_specs[[lipid]], shift)
}

individual_specs <- c(
  "dh_sph" = 0.65,
  "sph" = 1.20,
  "s1p" = 0.42,
  "dh_cer16_0" = 0.80,
  "dh_cer18_0" = 0.52,
  "dh_cer20_0" = 0.46,
  "dh_cer22_0" = 0.58,
  "dh_cer24_0" = 0.70,
  "dh_cer24_1" = 0.62,
  "cer16_0" = 1.45,
  "cer18_0" = 0.95,
  "cer20_0" = 0.84,
  "cer22_0" = 1.05,
  "cer24_0" = 1.30,
  "cer24_1" = 1.18,
  "dh_sm16_0" = 0.82,
  "dh_sm18_0" = 0.72,
  "dh_sm20_0" = 0.60,
  "dh_sm22_0" = 0.68,
  "dh_sm24_0" = 0.76,
  "dh_sm24_1" = 0.74,
  "sm16_0" = 5.20,
  "sm18_0" = 3.80,
  "sm20_0" = 2.90,
  "sm22_0" = 3.10,
  "sm24_0" = 3.40,
  "sm24_1" = 3.65,
  "hex_cer16_0" = 2.10,
  "hex_cer24_1" = 1.75,
  "lac_cer16_0" = 1.15,
  "lac_cer24_1" = 0.95
)

sms12_individual_demo <- design[c("infection", "condition")]
for (lipid in names(individual_specs)) {
  shift <- (match(lipid, names(individual_specs)) - 16) * 0.006
  sms12_individual_demo[[lipid]] <- make_measurement(individual_specs[[lipid]], shift)
}

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

make_pairwise_table <- function(lipids) {
  pairwise <- expand.grid(
    lipid = lipids,
    contrast = keep_contrasts,
    infection = infections,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  lipid_index <- match(pairwise$lipid, lipids)
  contrast_index <- match(pairwise$contrast, keep_contrasts)
  infection_direction <- ifelse(pairwise$infection == "yes", 1, -1)

  estimate <- (
    sin(lipid_index / 3) * 0.16 +
      cos(contrast_index / 2) * 0.12 +
      infection_direction * 0.08
  )

  p_value <- pmin(0.95, pmax(0.001, exp(-abs(estimate) * 12)))
  p_fdr <- pmin(0.99, p_value * 1.2)

  pairwise$estimate <- round(estimate, 4)
  pairwise$p.value <- signif(p_value, 3)
  pairwise$p_fdr <- signif(p_fdr, 3)
  pairwise
}

pairwise_total_demo <- make_pairwise_table(names(total_specs))
pairwise_individual_demo <- make_pairwise_table(names(individual_specs))

write_xlsx(
  list(
    sms12_total_demo = sms12_total_demo,
    sms12_individual_demo = sms12_individual_demo,
    pairwise_total_demo = pairwise_total_demo,
    pairwise_individual_demo = pairwise_individual_demo
  ),
  path = out_file
)

cat("Wrote demo workbook: ", out_file, "\n", sep = "")
cat("Rows in sms12_total_demo: ", nrow(sms12_total_demo), "\n", sep = "")
cat("Rows in sms12_individual_demo: ", nrow(sms12_individual_demo), "\n", sep = "")
cat("Rows in pairwise_total_demo: ", nrow(pairwise_total_demo), "\n", sep = "")
cat("Rows in pairwise_individual_demo: ", nrow(pairwise_individual_demo), "\n", sep = "")
