library(readxl)
library(janitor)
library(dplyr)
library(vegan)

# ==============================
# SMS1/2 individual lipid PERMANOVA and dispersion analysis
# - Reads one Excel sheet containing metadata and individual lipid columns
# - Uses PERMANOVA because individual lipid matrices can be too high-dimensional
#   for classical MANOVA relative to the number of samples
# - Applies log10(x + 1) transformation to lipid columns
# - Runs PERMANOVA for infection, condition, and their interaction
# - Runs betadisper/permutest to check group dispersion
# ==============================

## ---- Configure user-defined input and output paths ----
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"

## ---- Configure analysis options ----
metadata_cols <- c("infection", "condition")
infection_order <- c("no", "yes")
method <- "euclidean"
permutations <- 999
set_seed <- 1
log_transform <- TRUE

if (!file.exists(input_file)) {
  stop("Input file does not exist. Update input_file: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Read and prepare data ----
df_sms12_ind <- read_excel(input_file, sheet = sheet_name) |>
  clean_names()

missing_metadata <- setdiff(metadata_cols, names(df_sms12_ind))
if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata columns after clean_names(): ",
    paste(missing_metadata, collapse = ", ")
  )
}

lipid_cols <- setdiff(names(df_sms12_ind), metadata_cols)
if (length(lipid_cols) == 0) {
  stop("No lipid columns found after excluding metadata_cols.")
}

df_sms12_ind <- df_sms12_ind |>
  mutate(
    infection = factor(infection, levels = infection_order),
    condition = factor(condition)
  )

df_sms12_ind[lipid_cols] <- lapply(
  df_sms12_ind[lipid_cols],
  function(x) suppressWarnings(as.numeric(x))
)

if (log_transform) {
  has_negative <- any(
    vapply(df_sms12_ind[lipid_cols], function(x) any(x < 0, na.rm = TRUE), logical(1))
  )
  if (has_negative) {
    warning("Negative lipid values found. Check whether log10(x + 1) is appropriate.")
  }

  df_sms12_ind <- df_sms12_ind |>
    mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))
}

## ---- Lipid QC filtering ----
resp_mat <- as.matrix(df_sms12_ind[, lipid_cols, drop = FALSE])
storage.mode(resp_mat) <- "numeric"

all_na <- apply(resp_mat, 2, function(x) all(is.na(x)))
if (any(all_na)) {
  resp_mat <- resp_mat[, !all_na, drop = FALSE]
  lipid_cols <- colnames(resp_mat)
}

zero_var <- apply(resp_mat, 2, function(x) {
  x <- x[is.finite(x)]
  length(x) <= 1 || var(x) == 0
})
if (any(zero_var)) {
  resp_mat <- resp_mat[, !zero_var, drop = FALSE]
  lipid_cols <- colnames(resp_mat)
}

if (ncol(resp_mat) < 2) {
  stop("Fewer than two lipid columns remain after QC filtering.")
}

row_ok <- complete.cases(resp_mat, df_sms12_ind[, metadata_cols, drop = FALSE])
df_cc <- df_sms12_ind[row_ok, , drop = FALSE]
resp_mat_cc <- resp_mat[row_ok, , drop = FALSE]

df_cc <- df_cc |>
  mutate(
    infection = droplevels(infection),
    condition = droplevels(condition)
  )

if (nrow(df_cc) < 3) {
  stop("Fewer than three complete observations remain after filtering.")
}

if (nlevels(df_cc$infection) < 2) {
  stop("Fewer than two infection groups remain after filtering.")
}

if (nlevels(df_cc$condition) < 2) {
  stop("Fewer than two condition groups remain after filtering.")
}

cat(sprintf("Rows kept for analysis: %d / %d\n", nrow(df_cc), nrow(df_sms12_ind)))
cat(sprintf("Lipids kept after QC: %d\n", ncol(resp_mat_cc)))

writeLines(
  c(
    sprintf("Input file: %s", input_file),
    sprintf("Sheet name: %s", sheet_name),
    sprintf("Rows kept for analysis: %d / %d", nrow(df_cc), nrow(df_sms12_ind)),
    sprintf("Lipids kept after QC: %d", ncol(resp_mat_cc)),
    sprintf("Distance method: %s", method),
    sprintf("Permutations: %d", permutations),
    sprintf("Log10(x + 1) transformed: %s", log_transform)
  ),
  con = file.path(output_dir, "sms12_individual_lipids_analysis_summary.txt")
)

## ---- PERMANOVA ----
set.seed(set_seed)
sms12_ind_adonis <- adonis2(
  resp_mat_cc ~ infection * condition,
  data = df_cc,
  method = method,
  permutations = permutations,
  by = "term"
)

print(sms12_ind_adonis)

capture.output(
  sms12_ind_adonis,
  file = file.path(output_dir, "sms12_individual_lipids_permanova.txt")
)

## ---- Dispersion check ----
group_combined <- interaction(
  df_cc$infection,
  df_cc$condition,
  drop = TRUE
)

if (nlevels(group_combined) < 2) {
  stop("Fewer than two infection-condition groups remain after filtering.")
}

dist_ind <- dist(resp_mat_cc, method = method)

disp_ind <- betadisper(
  dist_ind,
  group_combined,
  type = "median",
  bias.adjust = TRUE
)

set.seed(set_seed)
disp_perm_ind <- permutest(disp_ind, permutations = permutations)

print(disp_perm_ind)

capture.output(
  disp_perm_ind,
  file = file.path(output_dir, "sms12_individual_lipids_dispersion_test.txt")
)
