library(readxl)
library(janitor)
library(dplyr)
library(tibble)
library(pheatmap)

# ==============================
# Average lipid heatmap script
# - Reads one Excel file
# - Cleans metadata and lipid columns
# - Applies log10(x + 1) transformation
# - Averages lipid values by cell line and infection
# - Z-scores each lipid across averaged groups
# - Creates a heatmap
# ==============================

## ---- Set your input and output here ----
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
out_dir    <- "path/to/your/output_directory/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Set factor order here ----
cell_line_order <- c("cell_line_1", "cell_line_2")
infection_order <- c("uninfected", "infected")

## ---- Read and clean data ----
df <- read_excel(input_file, sheet = sheet_name) %>%
  clean_names()

## ---- Define metadata and lipid columns ----
meta_cols  <- c("infection", "cell_line")
lipid_cols <- setdiff(names(df), meta_cols)

## ---- Force consistent factor ordering ----
df <- df %>%
  mutate(
    cell_line = factor(cell_line, levels = cell_line_order),
    infection = factor(infection, levels = infection_order)
  )

## ---- Ensure lipid columns are numeric ----
df[lipid_cols] <- lapply(df[lipid_cols], function(x) suppressWarnings(as.numeric(x)))

## ---- Log10 transform before averaging ----
df <- df %>%
  mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))

## ---- Average by cell line and infection ----
df_avg <- df %>%
  group_by(cell_line, infection) %>%
  summarise(
    across(all_of(lipid_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(sample_id = paste(cell_line, infection, sep = "_")) %>%
  arrange(cell_line, infection)

write.csv(
  df_avg,
  file.path(out_dir, "averaged_lipid_values.csv"),
  row.names = FALSE
)

## ---- Build matrix and z-score per lipid ----
mat <- df_avg %>%
  select(all_of(lipid_cols)) %>%
  as.matrix()

# Remove lipids with zero variance before scaling
keep_lipids <- apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)
mat <- mat[, keep_lipids, drop = FALSE]

mat <- scale(mat, center = TRUE, scale = TRUE)
rownames(mat) <- df_avg$sample_id

# pheatmap expects samples as columns here
mat_t <- t(mat)

## ---- Sample annotation ----
annotation_col <- df_avg %>%
  select(sample_id, infection, cell_line) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

stopifnot(identical(rownames(annotation_col), colnames(mat_t)))

## ---- Heatmap colour scale ----
pal <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(101)
bk  <- seq(-2, 2, length.out = length(pal) + 1)

## ---- Optional annotation colours ----
# Edit these if you want fixed annotation colours.
# The names must match the factor labels used in your data.
annotation_colours <- NULL

# Example:
# annotation_colours <- list(
#   infection = c("uninfected" = "#009E73", "infected" = "#CC79A7"),
#   cell_line = c("cell_line_1" = "#E69F00", "cell_line_2" = "#999999")
# )

## ---- Create heatmap ----
heatmap_obj <- pheatmap(
  mat_t,
  color = pal,
  breaks = bk,
  annotation_col = annotation_col,
  annotation_colors = annotation_colours,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = "white",
  filename = file.path(out_dir, "averaged_lipid_heatmap.png"),
  scale = "none",
  angle_col = 90,
  fontsize = 10
)

heatmap_obj
