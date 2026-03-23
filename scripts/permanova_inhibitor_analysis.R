
```r
library(readxl)
library(dplyr)
library(purrr)
library(janitor)
library(stringr)
library(tools)
library(vegan)

data_dir <- "path/to/your/input_directory/"   # Excel files (one per inhibitor)
METHOD   <- "euclidean"                       # "euclidean" for z-scored data
MIN_N    <- 4                                 # minimal observations for any subset analysis
out_dir  <- "path/to/your/output_directory/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read & align sheets (one file per inhibitor) ----
files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
stopifnot(length(files) >= 2)

read_inhibitor <- function(fp) {
  raw_name <- file_path_sans_ext(basename(fp))

  name_map <- c(
    "AKS466_new" = "AKS466",
    "DMSO_new" = "DMSO",
    "HPA-12_new" = "HPA-12",
    "Myriocin_new" = "Myriocin"
  )

  inhibitor_name <- unname(name_map[raw_name])

  read_excel(fp) |>
    janitor::clean_names() |>
    mutate(inhibitor = inhibitor_name)
}

dfs <- map(files, read_inhibitor)

# Union of columns to handle heterogeneity; keep consistent order
all_cols <- Reduce(union, lapply(dfs, names))
dfs_aligned <- lapply(dfs, function(dd) {
  missing <- setdiff(all_cols, names(dd))
  for (m in missing) dd[[m]] <- NA
  dd[, all_cols]
})

df <- bind_rows(dfs_aligned)

# ---- Factor cleaning ----
df <- df |>
  mutate(
    infection = factor(str_to_title(as.character(infection))),
    inhibitor = factor(as.character(inhibitor))
  )

# Optional: set readable orders (edit to your design)
if (all(c("Uninfected", "Infected") %in% levels(df$infection))) {
  df$infection <- factor(df$infection, levels = c("Uninfected", "Infected"))
}
df$inhibitor <- factor(df$inhibitor, levels = sort(unique(as.character(df$inhibitor))))

# ---- Identify lipid columns & coerce to numeric ----
meta_cols  <- c("infection", "inhibitor")
lipid_cols <- setdiff(names(df), meta_cols)

df <- df |>
  mutate(across(all_of(lipid_cols), ~ suppressWarnings(as.numeric(.))))

# ---- Log10 transform ----
df <- df |>
  mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))

# ---- Matrix + QC filtering ----
X <- as.matrix(df[, lipid_cols, drop = FALSE])

# Drop all-NA columns
all_na <- apply(X, 2, function(x) all(is.na(x)))
if (any(all_na)) {
  X <- X[, !all_na, drop = FALSE]
  lipid_cols <- colnames(X)
}

# Drop zero-variance columns (finite-only)
zero_var <- apply(X, 2, function(v) {
  v2 <- v[is.finite(v)]
  length(v2) <= 1 || sd(v2) == 0
})
if (any(zero_var)) {
  X <- X[, !zero_var, drop = FALSE]
  lipid_cols <- colnames(X)
}

# Scale (if Euclidean)
X_proc <- if (METHOD == "euclidean") scale(X) else X

# Keep complete rows
row_ok <- complete.cases(X_proc)
df_cc  <- df[row_ok, , drop = FALSE]
X_cc   <- X_proc[row_ok, , drop = FALSE]

df_cc <- df_cc %>%
  mutate(
    infection = factor(infection, levels = c("Uninfected", "Infected")),
    inhibitor = droplevels(inhibitor) # drop empty inhibitor levels if any
  )

cat(sprintf("Rows kept (complete cases): %d / %d\n", nrow(df_cc), nrow(df)))
cat(sprintf("Lipids kept after filtering: %d\n", ncol(X_cc)))

## =========================
## PERMANOVA (global, interaction, dispersion checks)
## =========================

# Global effects: PERMANOVA
# Analysis: PERMANOVA (adonis2) tested the effect of infection, inhibitor,
# and their interaction on the multivariate lipid profile.
# Permutation = 999, distance = Euclidean (on z-scores).
# Dispersion tests: betadisper + permutest checked whether group dispersions
# were homogeneous, to ensure PERMANOVA significance reflects centroid
# separation rather than spread.
#
# Output files:
# - inhibitors_ind_permanova_main.txt: infection and inhibitor effects
# - inhibitors_ind_permanova_dispersion_tests.txt: dispersion results for
#   infection and inhibitor groups

adon_main <- adonis2(
  X_cc ~ infection + inhibitor,
  data = df_cc,
  by = "term",
  method = METHOD,
  permutations = 999
)

print(adon_main)
capture.output(adon_main, file = file.path(out_dir, "inhibitors_ind_permanova_main.txt"))

# Dispersion checks (ensure significance isn't only spread differences)
D <- vegdist(X_cc, method = METHOD)

bd_inhib <- betadisper(D, df_cc$inhibitor)
bd_infec <- betadisper(D, df_cc$infection)

per_inhib <- permutest(bd_inhib)
per_infec <- permutest(bd_infec)

cat("\n-- Inhibitor betadisper permutest --\n")
print(per_inhib)

cat("\n-- Infection betadisper permutest --\n")
print(per_infec)
```
capture.output(
  list(inhibitor = per_inhib, infection = per_infec),
  file = file.path(out_dir, "inhibitors_ind_permanova_dispersion_tests.txt")
)