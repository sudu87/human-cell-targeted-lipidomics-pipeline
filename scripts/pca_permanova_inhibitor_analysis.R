```r
library(readxl)
library(janitor)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(tibble)
library(vegan)
library(tools)

# ==============================
# Clean PCA script
# - Reads one Excel file per inhibitor
# - Renames inhibitors immediately
# - Optionally subsets to Uninfected only
# - Performs log10 transform + QC + scaling + PCA + plot + PERMANOVA + dispersion checks
# ==============================

## ---- Set your directories here ----
data_dir <- "path/to/your/input_directory/"
METHOD   <- "euclidean"
out_dir  <- "path/to/your/output_directory/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Set whether to keep only Uninfected samples ----
KEEP_ONLY_UNINFECTED <- TRUE

## ---- Explicit filename -> inhibitor label mapping ----
name_map <- c(
  "AKS466_new" = "AKS466",
  "DMSO_new" = "DMSO",
  "HPA-12_new" = "HPA-12",
  "Myriocin_new" = "Myriocin"
)

## ---- Read all Excel files ----
files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
stopifnot(length(files) >= 1)

read_inhibitor <- function(fp, map_vec) {
  raw_name <- tools::file_path_sans_ext(basename(fp))

  if (!raw_name %in% names(map_vec)) {
    stop(
      "Filename not found in name_map: ", raw_name,
      "\nAdd it to name_map to keep naming reproducible."
    )
  }

  inhibitor_name <- unname(map_vec[raw_name])

  read_excel(fp) |>
    clean_names() |>
    mutate(inhibitor = inhibitor_name)
}

dfs <- purrr::map(files, read_inhibitor, map_vec = name_map)

## ---- Align columns across files and bind ----
all_cols <- Reduce(union, lapply(dfs, names))
dfs_aligned <- lapply(dfs, function(dd) {
  missing <- setdiff(all_cols, names(dd))
  for (m in missing) dd[[m]] <- NA
  dd[, all_cols]
})

df <- bind_rows(dfs_aligned)

## ---- Clean factors ----
df <- df |>
  mutate(
    infection = factor(str_to_title(as.character(infection))),
    inhibitor = factor(as.character(inhibitor), levels = unname(name_map))
  )

## ---- Optional: subset to Uninfected only ----
if (KEEP_ONLY_UNINFECTED) {
  df <- df |>
    filter(infection == "Uninfected") |>
    mutate(infection = droplevels(infection)) |>
    droplevels()
}

## ---- Identify lipid columns and coerce to numeric ----
meta_cols  <- c("infection", "inhibitor")
lipid_cols <- setdiff(names(df), meta_cols)

df <- df |>
  mutate(across(all_of(lipid_cols), ~ suppressWarnings(as.numeric(.))))

## ---- Log10 transform ----
df <- df |>
  mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))

## ---- Matrix and QC filtering ----
X <- as.matrix(df[, lipid_cols, drop = FALSE])

# Drop all-NA columns
all_na <- apply(X, 2, function(x) all(is.na(x)))
if (any(all_na)) {
  X <- X[, !all_na, drop = FALSE]
}

# Drop zero-variance columns
zero_var <- apply(X, 2, function(v) {
  v2 <- v[is.finite(v)]
  length(v2) <= 1 || sd(v2) == 0
})
if (any(zero_var)) {
  X <- X[, !zero_var, drop = FALSE]
}

lipid_cols <- colnames(X)

# Scale if Euclidean
X_proc <- if (METHOD == "euclidean") scale(X) else X

# Keep complete rows
row_ok <- complete.cases(X_proc)
df_cc  <- df[row_ok, , drop = FALSE]
X_cc   <- X_proc[row_ok, , drop = FALSE]

## ---- PCA ----
pca_fit <- prcomp(X_cc, center = TRUE, scale. = FALSE)

# Variance explained
pca_var <- pca_fit$sdev^2
pca_var_exp <- pca_var / sum(pca_var)
pca_var_df <- tibble(
  PC = paste0("PC", seq_along(pca_var)),
  Variance = pca_var,
  Pct = 100 * pca_var_exp
)

write.csv(
  pca_var_df,
  file.path(out_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

# Scores
scores <- as.data.frame(pca_fit$x) |>
  mutate(
    infection = df_cc$infection,
    inhibitor = df_cc$inhibitor
  )

var_exp <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))

## ---- PCA plot ----
p_pca_2 <- ggplot(scores, aes(x = PC1, y = PC2, fill = inhibitor)) +
  geom_point(size = 3, alpha = 0.9, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c(
    "AKS466" = "#5D3A9B",
    "DMSO" = "#E66100",
    "HPA-12" = "#0C7BDC",
    "Myriocin" = "#FFC20A"
  )) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2])
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, colour = "black"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

## ---- PERMANOVA ----
adon_uninfected_ind_dmso <- adonis2(
  X_cc ~ inhibitor,
  data = df_cc,
  by = "term",
  method = METHOD,
  permutations = 999
)

print(adon_uninfected_ind_dmso)

capture.output(
  adon_uninfected_ind_dmso,
  file = file.path(out_dir, "inhibitors_total_permanova_main.txt")
)

## ---- Dispersion checks ----
D <- vegdist(X_cc, method = METHOD)

bd_inhib <- betadisper(D, df_cc$inhibitor)
per_inhib <- permutest(bd_inhib)

cat("\n-- Inhibitor betadisper permutest --\n")
print(per_inhib)

capture.output(
  list(inhibitor = per_inhib),
  file = file.path(out_dir, "inhibitors_total_permanova_dispersion_tests.txt")
)
````
