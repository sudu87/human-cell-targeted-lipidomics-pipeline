library(readxl)
library(janitor)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
library(vegan)

# ==============================
# PCA and PERMANOVA Hela HUVEC analysis
# - Reads one Excel file
# - Uses infection and cell_line as metadata columns
# - Performs log10 transform + PERMANOVA + PCA + PCA plots
# ==============================

## ---- Set your file paths here ----
input_file <- "path/to/your/input_file.xlsx"
out_dir    <- "path/to/your/output_directory/"
METHOD     <- "euclidean"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Read the Excel file ----
df <- read_excel(input_file) |>
  clean_names() |>
  mutate(
    cell_line = factor(str_to_title(as.character(cell_line))),
    infection = factor(str_to_title(as.character(infection)))
  )

## ---- Define lipid columns ----
meta_cols  <- c("infection", "cell_line")
lipid_cols <- setdiff(names(df), meta_cols)

## ---- Log10 transform lipid columns ----
df <- df |>
  mutate(across(all_of(lipid_cols), ~ log10(as.numeric(.x) + 1)))

## ---- Build lipid matrix ----
idx <- complete.cases(df[, lipid_cols], df$infection, df$cell_line)

X <- as.matrix(df[idx, lipid_cols, drop = FALSE])

meta <- df[idx, c("cell_line", "infection")] |>
  droplevels()

## ---- Remove zero-variance lipids ----
nzv <- vapply(as.data.frame(X), function(x) var(x, na.rm = TRUE) > 0, logical(1))

X <- X[, nzv, drop = FALSE]
lipid_cols <- colnames(X)

## ---- PERMANOVA ----
set.seed(1)

permanova_fit <- adonis2(
  X ~ infection * cell_line,
  data = meta,
  method = METHOD,
  permutations = 999,
  by = "term"
)

print(permanova_fit)

capture.output(
  permanova_fit,
  file = file.path(out_dir, "permanova_results.txt")
)

## ---- PCA ----
pca_fit <- prcomp(X, center = TRUE, scale. = TRUE)

pca_var <- pca_fit$sdev^2
pca_var_exp <- pca_var / sum(pca_var)

pca_var_df <- tibble(
  PC  = paste0("PC", seq_along(pca_var)),
  Var = pca_var,
  Pct = 100 * pca_var_exp
)

write.csv(
  pca_var_df,
  file.path(out_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

scores <- as_tibble(pca_fit$x, .name_repair = "minimal") |>
  bind_cols(meta)

var_exp <- 100 * pca_var_exp

## ---- PCA plot ----
p_pca_2 <- ggplot(
  scores,
  aes(x = PC1, y = PC2, fill = infection, shape = cell_line)
) +
  geom_point(
    size = 3,
    alpha = 0.9,
    color = "black",
    stroke = 0.8
  ) +
  scale_shape_manual(values = c(21, 24, 22, 23, 25)) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2])
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, colour = "black")
    ),
    shape = guide_legend(
      override.aes = list(fill = "white", colour = "black")
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

ggsave(
  file.path(out_dir, "pca_pc1_pc2.png"),
  p_pca_2,
  width = 7,
  height = 5,
  dpi = 300
)

## ---- Scree plot ----
p_scree <- ggplot(pca_var_df, aes(x = seq_along(Pct), y = Pct)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = "PCA Scree Plot"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "pca_scree_plot.png"),
  p_scree,
  width = 7,
  height = 4,
  dpi = 300
)
