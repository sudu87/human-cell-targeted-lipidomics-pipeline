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
# PCA and PERMANOVA analysis for total sphingolipids
# - Reads one Excel file per inhibitor/control condition
# - Keeps water-soluble and DMSO-soluble analyses separate
# - Applies log10(x + 1), QC filtering, z-scoring, PCA, PERMANOVA,
#   and dispersion checks
# ==============================

## ---- Choose analysis group ----
# Use "dmso" for DMSO-soluble inhibitors and DMSO control.
# Use "untreated" for water-soluble inhibitors and untreated control.
ANALYSIS_GROUP <- "dmso"

METHOD <- "euclidean"
MIN_N  <- 4

## ---- Configure input files and labels ----
analysis_configs <- list(
  dmso = list(
    data_dir = "path/to/your/dmso_input_directory",
    out_dir = "path/to/your/dmso_output_directory",
    name_map = c(
      "control_file_name_without_extension" = "Control label",
      "condition_file_name_without_extension" = "Condition label"
    )
  ),
  untreated = list(
    data_dir = "path/to/your/untreated_input_directory",
    out_dir = "path/to/your/untreated_output_directory",
    name_map = c(
      "control_file_name_without_extension" = "Control label",
      "condition_file_name_without_extension" = "Condition label"
    )
  )
)

if (!ANALYSIS_GROUP %in% names(analysis_configs)) {
  stop(
    "Unknown ANALYSIS_GROUP: ", ANALYSIS_GROUP,
    "\nUse one of: ", paste(names(analysis_configs), collapse = ", ")
  )
}

cfg <- analysis_configs[[ANALYSIS_GROUP]]
data_dir <- cfg$data_dir
out_dir  <- cfg$out_dir
name_map <- cfg$name_map

if (!dir.exists(data_dir)) {
  stop(
    "Input directory does not exist. Update data_dir for ANALYSIS_GROUP '",
    ANALYSIS_GROUP,
    "': ",
    data_dir
  )
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Read and align sheets: one file per inhibitor/control ----
files <- list.files(data_dir, pattern = "\\.xlsx$", full.names = TRUE)
if (length(files) < 2) {
  stop("Expected at least two .xlsx files in data_dir: ", data_dir)
}

read_inhibitor <- function(fp, map_vec) {
  raw_name <- file_path_sans_ext(basename(fp))

  if (!raw_name %in% names(map_vec)) {
    stop(
      "Filename not found in name_map: ", raw_name,
      "\nAdd this file name to name_map before running the analysis."
    )
  }

  inhibitor_name <- unname(map_vec[raw_name])

  read_excel(fp) |>
    janitor::clean_names() |>
    mutate(inhibitor = inhibitor_name)
}

dfs <- map(files, read_inhibitor, map_vec = name_map)

# Union of columns handles lipids that are present in some files but absent in others.
all_cols <- Reduce(union, lapply(dfs, names))
dfs_aligned <- lapply(dfs, function(dd) {
  missing <- setdiff(all_cols, names(dd))
  for (m in missing) dd[[m]] <- NA
  dd[, all_cols]
})

df <- bind_rows(dfs_aligned)

## ---- Factor cleaning ----
df <- df |>
  mutate(
    infection = factor(str_to_title(as.character(infection))),
    inhibitor = factor(as.character(inhibitor), levels = unname(name_map))
  )

if (all(c("Uninfected", "Infected") %in% levels(df$infection))) {
  df$infection <- factor(df$infection, levels = c("Uninfected", "Infected"))
}

df$inhibitor <- droplevels(df$inhibitor)

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

# Drop all-NA lipid columns.
all_na <- apply(X, 2, function(x) all(is.na(x)))
if (any(all_na)) {
  X <- X[, !all_na, drop = FALSE]
  lipid_cols <- colnames(X)
}

# Drop lipid columns with zero variance among finite values.
zero_var <- apply(X, 2, function(v) {
  v2 <- v[is.finite(v)]
  length(v2) <= 1 || sd(v2) == 0
})
if (any(zero_var)) {
  X <- X[, !zero_var, drop = FALSE]
  lipid_cols <- colnames(X)
}

if (ncol(X) < 2) {
  stop("Fewer than two lipid columns remain after QC filtering.")
}

# Scale if Euclidean distance is used.
X_proc <- if (METHOD == "euclidean") scale(X) else X

# Keep complete rows.
row_ok <- complete.cases(X_proc)
df_cc  <- df[row_ok, , drop = FALSE]
X_cc   <- X_proc[row_ok, , drop = FALSE]

df_cc <- df_cc |>
  mutate(
    infection = droplevels(infection),
    inhibitor = droplevels(inhibitor)
  )

cat(sprintf("Analysis group: %s\n", ANALYSIS_GROUP))
cat(sprintf("Rows kept (complete cases): %d / %d\n", nrow(df_cc), nrow(df)))
cat(sprintf("Lipids kept after filtering: %d\n", ncol(X_cc)))

if (nrow(df_cc) < MIN_N) {
  stop("Fewer than MIN_N complete observations remain after filtering.")
}

if (nlevels(df_cc$inhibitor) < 2) {
  stop("Fewer than two inhibitor/control groups remain after filtering.")
}

## ---- PCA ----
pca_fit <- prcomp(X_cc, center = TRUE, scale. = FALSE)

pca_var <- pca_fit$sdev^2
pca_var_exp <- pca_var / sum(pca_var)
pca_var_df <- tibble(
  PC = paste0("PC", seq_along(pca_var)),
  Variance = pca_var,
  Pct = 100 * pca_var_exp
)

write.csv(
  pca_var_df,
  file.path(out_dir, paste0(ANALYSIS_GROUP, "_pca_variance_explained.csv")),
  row.names = FALSE
)

scores <- as.data.frame(pca_fit$x) |>
  mutate(
    infection = df_cc$infection,
    inhibitor = df_cc$inhibitor
  )

var_exp <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))

## ---- PCA plot ----
p_pca <- ggplot(scores, aes(x = PC1, y = PC2, fill = inhibitor)) +
  geom_point(size = 3, alpha = 0.9, shape = 21, color = "black", stroke = 0.8) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    fill = "Condition"
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

ggsave(
  filename = file.path(out_dir, paste0(ANALYSIS_GROUP, "_pca_plot.pdf")),
  plot = p_pca,
  width = 6.5,
  height = 5,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, paste0(ANALYSIS_GROUP, "_pca_plot.png")),
  plot = p_pca,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 300
)

## ---- PERMANOVA ----
adon_inhibitor <- adonis2(
  X_cc ~ inhibitor,
  data = df_cc,
  by = "term",
  method = METHOD,
  permutations = 999
)

print(adon_inhibitor)

capture.output(
  adon_inhibitor,
  file = file.path(out_dir, paste0(ANALYSIS_GROUP, "_permanova_main.txt"))
)

## ---- Dispersion checks ----
D <- vegdist(X_cc, method = METHOD)

bd_inhib <- betadisper(D, df_cc$inhibitor)
per_inhib <- permutest(bd_inhib)

cat("\n-- Inhibitor betadisper permutest --\n")
print(per_inhib)

capture.output(
  list(inhibitor = per_inhib),
  file = file.path(out_dir, paste0(ANALYSIS_GROUP, "_permanova_dispersion_tests.txt"))
)
