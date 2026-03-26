
```r
## ---- Set your output directory here if needed ----
out_dir <- "path/to/your/output_directory/"

## ---- Create diagnostics output folder ----
out_dir_diag <- file.path(out_dir, "diagnostics")
dir.create(out_dir_diag, showWarnings = FALSE, recursive = TRUE)

## ---- Extract all lipid values into one vector ----
all_values <- unlist(df[lipid_cols])

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
  log10(all_values + 1),
  breaks = 40,
  col = "tomato",
  border = "white",
  main = "Distribution after log10 transformation",
  xlab = "log10(Lipid concentration + 1)"
)
dev.off()

## Density plot of raw values
png(file.path(out_dir_diag, "density_raw_values.png"), width = 1000, height = 800, res = 150)
plot(
  density(all_values, na.rm = TRUE),
  main = "Density (raw)",
  xlab = "Raw values",
  col = "steelblue",
  lwd = 2
)
dev.off()

## Density plot of log10-transformed values
png(file.path(out_dir_diag, "density_log10_values.png"), width = 1000, height = 800, res = 150)
plot(
  density(log10(all_values + 1), na.rm = TRUE),
  main = "Density (log10)",
  xlab = "log10(values + 1)",
  col = "tomato",
  lwd = 2
)
dev.off()
## ---- QQ-plot of pooled ANOVA residuals ----
all_resids <- unlist(lapply(lipid_cols, function(lip) {
  form <- reformulate("infection", response = lip)
  fit <- aov(form, data = df)
  residuals(fit)
}))

png(file.path(out_dir_diag, "qqplot_residuals.png"), width = 1000, height = 800, res = 150)
qqnorm(all_resids, main = "QQ-plot of pooled ANOVA residuals")
qqline(all_resids, col = "red", lwd = 2)
dev.off()
