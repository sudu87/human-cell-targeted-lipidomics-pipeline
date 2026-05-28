library(readxl)
library(janitor)
library(dplyr)
library(pheatmap)

# ==============================
# SMS1/2 ordered-condition lipid heatmap
# - Reads one raw Excel sheet containing infection, condition, and lipid columns
# - Subsets by infection status
# - Applies log10(x + 1) transformation to the selected subset
# - Averages replicates within each condition
# - Draws row-scaled lipid heatmaps using a user-defined condition order
# ==============================

## ---- Configure user-defined input and output paths ----
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- "your_sheet_name"
output_dir <- "path/to/your/output_directory"

## ---- Configure subset and ordering ----
# Use c("yes") for infected only, c("no") for uninfected only,
# or c("yes", "no") to create one heatmap for each subset.
infection_groups_to_plot <- c("yes", "no")

condition_order <- c(
  "sms12_WT",
  "sms12_DKO",
  "sms1_16_no_AHT",
  "sms1_16_AHT",
  "sms2_no_AHT",
  "sms2_AHT",
  "sms2-M64R_no_AHT",
  "sms2-M64R_AHT"
)

condition_labels <- c(
  "sms12_WT" = "HeLa WT (-AHT)",
  "sms12_DKO" = "SMS1/2 DKO (-AHT)",
  "sms1_16_no_AHT" = "SMS1 comp (-AHT)",
  "sms1_16_AHT" = "SMS1 comp (+AHT)",
  "sms2_no_AHT" = "SMS2 comp (-AHT)",
  "sms2_AHT" = "SMS2 comp (+AHT)",
  "sms2-M64R_no_AHT" = "SMS2-M64R (-AHT)",
  "sms2-M64R_AHT" = "SMS2-M64R (+AHT)"
)

output_prefix <- "sms12_ordered_condition_heatmap"
log_transform <- TRUE
cluster_rows <- FALSE
cluster_cols <- FALSE

## ---- Read and prepare raw data ----
if (!file.exists(input_file)) {
  stop("Input file does not exist. Update input_file: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_excel(input_file, sheet = sheet_name) |>
  clean_names()

meta_cols <- c("infection", "condition")
missing_meta <- setdiff(meta_cols, names(df))
if (length(missing_meta) > 0) {
  stop(
    "Missing required columns after clean_names(): ",
    paste(missing_meta, collapse = ", ")
  )
}

lipid_cols <- setdiff(names(df), meta_cols)
if (length(lipid_cols) == 0) {
  stop("No lipid columns found after excluding metadata columns.")
}

## ---- Heatmap helper ----
make_heatmap_for_infection <- function(infection_value) {
  df_subset <- df |>
    filter(infection == infection_value)

  if (nrow(df_subset) == 0) {
    warning("No rows found for infection value: ", infection_value)
    return(NULL)
  }

  df_subset[lipid_cols] <- lapply(
    df_subset[lipid_cols],
    function(x) suppressWarnings(as.numeric(x))
  )

  if (log_transform) {
    has_negative <- any(
      vapply(df_subset[lipid_cols], function(x) any(x < 0, na.rm = TRUE), logical(1))
    )
    if (has_negative) {
      warning(
        "Negative lipid values found for infection value '",
        infection_value,
        "'. Check whether log10(x + 1) is appropriate."
      )
    }

    df_subset <- df_subset |>
      mutate(across(all_of(lipid_cols), ~ log10(.x + 1)))
  }

  avg <- df_subset |>
    filter(!is.na(condition)) |>
    group_by(condition) |>
    summarise(
      across(all_of(lipid_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(across(all_of(lipid_cols), ~ ifelse(is.nan(.x), NA_real_, .x)))

  if (nrow(avg) == 0) {
    warning("No non-missing conditions found for infection value: ", infection_value)
    return(NULL)
  }

  infection_file_label <- gsub("[^A-Za-z0-9_+-]+", "_", infection_value)

  write.csv(
    avg,
    file = file.path(
      output_dir,
      paste0(output_prefix, "_", infection_file_label, "_averaged_values.csv")
    ),
    row.names = FALSE
  )

  mat <- t(as.matrix(avg[, lipid_cols, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  colnames(mat) <- as.character(avg$condition)

  keep <- condition_order[condition_order %in% colnames(mat)]
  extras <- setdiff(colnames(mat), keep)
  mat <- mat[, c(keep, extras), drop = FALSE]

  row_has_values <- apply(mat, 1, function(x) any(is.finite(x)))
  mat <- mat[row_has_values, , drop = FALSE]

  if (nrow(mat) == 0) {
    warning("No lipids with finite values for infection value: ", infection_value)
    return(NULL)
  }

  if (ncol(mat) > 1) {
    row_var <- apply(mat, 1, function(x) {
      x <- x[is.finite(x)]
      length(x) > 1 && var(x) > 0
    })
    mat <- mat[row_var, , drop = FALSE]
  }

  if (nrow(mat) == 0) {
    warning(
      "No lipids with non-zero variance for infection value: ",
      infection_value
    )
    return(NULL)
  }

  display_labels <- condition_labels[colnames(mat)]
  missing_labels <- is.na(display_labels)
  display_labels[missing_labels] <- colnames(mat)[missing_labels]
  colnames(mat) <- display_labels

  pal <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(101)
  breaks <- seq(-2, 2, length.out = length(pal) + 1)

  scale_mode <- if (ncol(mat) > 1) "row" else "none"

  heatmap_args <- list(
    mat = mat,
    scale = scale_mode,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    na_col = "grey90",
    color = pal,
    border_color = "white",
    fontsize = 10,
    angle_col = 90,
    breaks = breaks,
    silent = TRUE
  )

  heatmap_png <- do.call(
    pheatmap,
    c(
      heatmap_args,
      list(
        filename = file.path(
          output_dir,
          paste0(output_prefix, "_", infection_file_label, ".png")
        ),
        width = 8,
        height = 8
      )
    )
  )

  heatmap_pdf <- do.call(
    pheatmap,
    c(
      heatmap_args,
      list(
        filename = file.path(
          output_dir,
          paste0(output_prefix, "_", infection_file_label, ".pdf")
        ),
        width = 8,
        height = 8
      )
    )
  )

  list(
    infection = infection_value,
    averaged_values = avg,
    matrix = mat,
    heatmap_png = heatmap_png,
    heatmap_pdf = heatmap_pdf
  )
}

heatmap_results <- lapply(infection_groups_to_plot, make_heatmap_for_infection)
names(heatmap_results) <- infection_groups_to_plot
