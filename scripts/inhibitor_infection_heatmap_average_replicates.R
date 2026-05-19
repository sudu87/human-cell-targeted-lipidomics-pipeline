library(readxl)
library(janitor)
library(dplyr)
library(tibble)
library(pheatmap)

if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

## ---- Set your file here each time ----
file_path <- "path/to/your/input_file.xlsx"

## ---- Load & clean ----
df <- read_excel(file_path) %>%
  clean_names()

meta_cols <- c("infection", "inhibitor")
lipid_cols <- setdiff(names(df), meta_cols)

## ---- Factor ordering (consistent) ----
df <- df %>%
  mutate(
    inhibitor = factor(
      inhibitor,
      levels = c("untreated", "DMSO", "AKS466", "HPA-12", "Desipramin", "Myriocin", "ARC39")
    ),
    infection = factor(
      infection,
      levels = c("Uninfected", "infected")
    )
  )

## ---- Make sure lipid columns are numeric ----
df[lipid_cols] <- lapply(df[lipid_cols], function(x) suppressWarnings(as.numeric(x)))

## ---- Log10 transform BEFORE averaging ----
df <- df %>%
  mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))

## ---- Average replicates & enforce sample order ----
df_avg <- df %>%
  group_by(inhibitor, infection) %>%
  summarise(
    across(all_of(lipid_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(sample_id = paste(inhibitor, infection, sep = "_")) %>%
  mutate(
    order_key = case_when(
      infection == "Uninfected" & inhibitor %in% c("untreated", "DMSO") ~ 1,
      infection == "Uninfected" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(order_key, inhibitor, infection) %>%
  dplyr::select(-order_key)

## ---- Build matrix & z-score ----
mat <- df_avg %>%
  dplyr::select(all_of(lipid_cols)) %>%
  as.matrix()

mat <- scale(mat, center = TRUE, scale = TRUE)
rownames(mat) <- df_avg$sample_id

mat_t <- t(mat)

annotation_col <- df_avg %>%
  dplyr::select(sample_id, infection, inhibitor) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

stopifnot(identical(rownames(annotation_col), colnames(mat_t)))

## ---- Heatmap palettes ----
pal <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(101)

ann_colors <- list(
  infection = c(
    "Uninfected" = "#009E73",
    "infected" = "#CC79A7"
  )
)

breaks <- seq(-2, 2, length.out = 101)

## ---- Make heatmap object ----
heatmap_object <- pheatmap(
  mat_t,
  color = pal,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = "white",
  filename = NA,
  scale = "none",
  angle_col = 90,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  breaks = breaks
)
