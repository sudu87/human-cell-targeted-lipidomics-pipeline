library(readxl)
library(janitor)

## ---- Configure user-defined input and output paths ----
input_file <- "path/to/your/input_file.xlsx"
sheet_name <- NULL
out_dir <- "path/to/your/output_directory"

# Add any non-lipid metadata columns in your raw data to this vector.
metadata_cols <- c(
  "infection",
  "condition",
  "inhibitor",
  "cell_line",
  "sample",
  "sample_id",
  "sample_name",
  "replicate",
  "group"
)
infection_col <- "infection"

if (!file.exists(input_file)) {
  stop("Input file does not exist. Update input_file: ", input_file)
}

## ---- Read and prepare raw data ----
df <- if (is.null(sheet_name) || identical(sheet_name, "")) {
  read_excel(input_file)
} else {
  read_excel(input_file, sheet = sheet_name)
}

df <- clean_names(df)

if (!infection_col %in% names(df)) {
  stop("The cleaned input data must contain an infection column.")
}

metadata_cols <- unique(c(infection_col, metadata_cols))
candidate_lipid_cols <- setdiff(names(df), metadata_cols)

if (length(candidate_lipid_cols) == 0) {
  stop("No candidate lipid columns found after excluding metadata_cols.")
}

df[candidate_lipid_cols] <- lapply(
  df[candidate_lipid_cols],
  function(x) suppressWarnings(as.numeric(x))
)

lipid_cols <- candidate_lipid_cols[
  vapply(df[candidate_lipid_cols], function(x) any(is.finite(x)), logical(1))
]

if (length(lipid_cols) == 0) {
  stop("No numeric lipid columns found. Check metadata_cols and the input data.")
}

dropped_cols <- setdiff(candidate_lipid_cols, lipid_cols)
if (length(dropped_cols) > 0) {
  message(
    "Dropped non-numeric or empty columns: ",
    paste(dropped_cols, collapse = ", ")
  )
}

df[[infection_col]] <- factor(df[[infection_col]])

cat(sprintf("Rows loaded: %d\n", nrow(df)))
cat(sprintf("Lipids included: %d\n", length(lipid_cols)))

## ---- Create diagnostics output folder ----
out_dir_diag <- file.path(out_dir, "diagnostics")
dir.create(out_dir_diag, showWarnings = FALSE, recursive = TRUE)

## ---- Extract all lipid values into one vector ----
all_values <- unlist(df[lipid_cols], use.names = FALSE)
all_values <- all_values[is.finite(all_values)]

if (length(all_values) == 0) {
  stop("No finite lipid values found for diagnostics.")
}

if (any(all_values < 0)) {
  warning("Negative lipid values found. Check whether log10(x + 1) is appropriate.")
}

log_values <- log10(all_values + 1)
log_values <- log_values[is.finite(log_values)]

if (length(log_values) == 0) {
  stop("No finite log10-transformed lipid values found for diagnostics.")
}

## ---- Histogram of raw values ----
png(file.path(out_dir_diag, "hist_raw_values.png"), width = 1000, height = 800, res = 150)
hist(
  all_values,
  breaks = 40,
  col = "steelblue",
  border = "white",
  main = "Distribution of raw lipid concentrations",
  xlab = "Lipid concentration"
)
dev.off()

## ---- Histogram of log10-transformed values ----
png(file.path(out_dir_diag, "hist_log10_values.png"), width = 1000, height = 800, res = 150)
hist(
  log_values,
  breaks = 40,
  col = "tomato",
  border = "white",
  main = "Distribution after log10 transformation",
  xlab = "log10(Lipid concentration + 1)"
)
dev.off()

plot_density <- function(values, file, main, xlab, col) {
  png(file, width = 1000, height = 800, res = 150)
  if (length(unique(values)) < 2) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, "Density plot requires at least two distinct values.")
  } else {
    plot(
      density(values),
      main = main,
      xlab = xlab,
      col = col,
      lwd = 2
    )
  }
  dev.off()
}

## Density plot of raw values
plot_density(
  all_values,
  file.path(out_dir_diag, "density_raw_values.png"),
  "Density (raw)",
  "Raw values",
  "steelblue"
)

## Density plot of log10-transformed values
plot_density(
  log_values,
  file.path(out_dir_diag, "density_log10_values.png"),
  "Density (log10)",
  "log10(values + 1)",
  "tomato"
)

## ---- QQ-plot of pooled ANOVA residuals ----
all_resids <- unlist(lapply(lipid_cols, function(lip) {
  model_df <- df[, c(infection_col, lip), drop = FALSE]
  model_df <- model_df[complete.cases(model_df), , drop = FALSE]
  model_df[[infection_col]] <- droplevels(factor(model_df[[infection_col]]))

  if (nrow(model_df) < 3 || nlevels(model_df[[infection_col]]) < 2) {
    return(numeric(0))
  }

  form <- reformulate(infection_col, response = lip)
  fit <- aov(form, data = model_df)
  residuals(fit)
}), use.names = FALSE)

if (length(all_resids) >= 2) {
  png(file.path(out_dir_diag, "qqplot_residuals.png"), width = 1000, height = 800, res = 150)
  qqnorm(all_resids, main = "QQ-plot of pooled ANOVA residuals")
  qqline(all_resids, col = "red", lwd = 2)
  dev.off()
} else {
  warning("Skipping QQ plot because too few ANOVA residuals could be calculated.")
}
