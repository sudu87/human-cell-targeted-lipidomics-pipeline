library(readxl)
library(janitor)
library(dplyr)
library(ggplot2)
library(tibble)
library(vegan)

# ==============================
# SMS1/2 total lipid MANOVA, PCA, and dispersion checks
# - Reads one Excel sheet containing metadata and lipid abundance columns
# - Applies log10(x + 1) transformation to lipid columns
# - Runs MANOVA, PCA, scree/PCA plots, and betadisper dispersion checks
# ==============================

## ---- Configure input and outputs ----
input_file <- "sms12_lipidomics/SMS1_2_Total lipids_no_19_test_file.xlsx"
sheet_name <- "SMS1&2_Total lipids_uninfected"

out_dir <- "sms12_lipidomics/statistics"
plot_dir <- "sms12_lipidomics"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

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

## ---- Matrix and complete-case filtering ----
resp_mat <- as.matrix(sms_test[, lipid_cols, drop = FALSE])

row_ok <- complete.cases(resp_mat, sms_test[, meta_cols, drop = FALSE])
sms_test_cc <- sms_test[row_ok, , drop = FALSE]
resp_mat_cc <- resp_mat[row_ok, , drop = FALSE]

sms_test_cc <- sms_test_cc |>
  mutate(
    infection = droplevels(infection),
    condition = droplevels(condition)
  )

cat(sprintf("Rows kept for analysis: %d / %d\n", nrow(sms_test_cc), nrow(sms_test)))
cat(sprintf("Lipids included: %d\n", ncol(resp_mat_cc)))

## ---- MANOVA ----
sms_test_manova <- manova(resp_mat_cc ~ infection * condition, data = sms_test_cc)
manova_pillai <- summary(sms_test_manova, test = "Pillai")

print(manova_pillai)

capture.output(
  manova_pillai,
  file = file.path(out_dir, "manova_pillai.txt")
)

## ---- Colors ----
cond_levels_test <- levels(droplevels(sms_test_cc$condition))

okabe_ito <- c(
  "#0072B2", # blue
  "#56B4E9", # sky blue
  "#D73027", # red
  "#CC79A7", # reddish purple
  "#D55E00", # vermillion
  "#F6C667", # light orange
  "#006D5B", # dark green
  "#66C2A5"  # light green
)

condition_cols_test <- setNames(
  okabe_ito[seq_along(cond_levels_test)],
  cond_levels_test
)

## ---- PCA ----
pca_sms_test <- prcomp(resp_mat_cc, center = TRUE, scale. = TRUE)

pca_var_test <- pca_sms_test$sdev^2
pca_var_exp_test <- pca_var_test / sum(pca_var_test)

pca_var_test_df <- tibble(
  PC = paste0("PC", seq_along(pca_var_test)),
  Var = pca_var_test,
  Pct = 100 * pca_var_exp_test
)

write.csv(
  pca_var_test_df,
  file = file.path(plot_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

pve_test <- summary(pca_sms_test)$importance[2, 1:2] * 100
pc1_lab_test <- paste0("PC1 (", round(pve_test[1], 1), "%)")
pc2_lab_test <- paste0("PC2 (", round(pve_test[2], 1), "%)")

scores_sms_test <- as.data.frame(pca_sms_test$x[, 1:2, drop = FALSE]) |>
  mutate(
    infection = sms_test_cc$infection,
    condition = sms_test_cc$condition
  )

write.csv(
  scores_sms_test,
  file = file.path(plot_dir, "pca_scores.csv"),
  row.names = FALSE
)

## ---- PCA scree plot ----
p_scree_sms12_test <- ggplot(pca_var_test_df, aes(x = seq_along(Pct), y = Pct)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = "PCA Scree Plot"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(plot_dir, "pca_scree.pdf"),
  plot = p_scree_sms12_test,
  width = 7,
  height = 4,
  dpi = 300
)

## ---- PCA score plot ----
p_pca_sms12_test <- ggplot(
  scores_sms_test,
  aes(x = PC1, y = PC2, color = condition, shape = infection)
) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = condition_cols_test) +
  labs(
    x = pc1_lab_test,
    y = pc2_lab_test,
    color = "Condition",
    shape = "Infection",
    title = "PCA of SMS1/2 total lipids"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

ggsave(
  filename = file.path(plot_dir, "pca_scores.pdf"),
  plot = p_pca_sms12_test,
  width = 6.5,
  height = 5,
  dpi = 300
)

## ---- Dispersion checks ----
group_combined_test <- interaction(
  sms_test_cc$infection,
  sms_test_cc$condition,
  drop = TRUE
)

D_test <- dist(resp_mat_cc, method = "euclidean")

disp_groups_test <- betadisper(
  D_test,
  group_combined_test,
  type = "median",
  bias.adjust = TRUE
)

set.seed(1)
per_groups_test <- permutest(disp_groups_test, permutations = 999)

print(per_groups_test)

capture.output(
  per_groups_test,
  file = file.path(out_dir, "dispersion_test.txt")
)
