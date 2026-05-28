library(readxl)
library(dplyr)
library(ggplot2)

## ---- Set your file here each time ----
file_path <- "path/to/your/input_file.xlsx"

## ---- Load data ----
df_heat <- read_excel(file_path)

## ---- Inspect available contrasts ----
dput(unique(df_heat$contrast))

## ---- Clean and transform ----
df_hm <- df_heat %>%
  filter(
    !is.na(lipid),
    !is.na(estimate),
    !is.na(p.value),
    !is.na(infection),
    !is.na(contrast)
  ) %>%
  mutate(
    log2FC = log2(10^estimate),
    log2FC_sig = ifelse(p.value < 0.05, log2FC, NA_real_)
  )

## ---- Set custom lipid order here ----
lipid_order <- c(
  "lipid1",
  "lipid2"
)

## ---- Contrast ordering ----
contrast_levels <- df_hm %>%
  distinct(contrast) %>%
  arrange(contrast) %>%
  pull(contrast)

## ---- Apply ordering ----
plot_df <- df_hm %>%
  mutate(
    lipid = factor(lipid, levels = rev(lipid_order)),
    contrast = factor(contrast, levels = contrast_levels),
    infection = factor(infection, levels = c("Uninfected", "infected"))
  )

## ---- Set legend title here ----
legend_title <- "Treatment\nvs\nUntreated"

## ---- Heatmap ----
heatmap_plot <- ggplot(plot_df, aes(x = contrast, y = lipid, fill = log2FC_sig)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_y_discrete(position = "right") +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey85",
    name = legend_title,
    guide = guide_colorbar(barwidth = 0.7)
  ) +
  facet_grid(~ infection) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 0),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10, angle = 90),
    legend.title = element_text(size = 10, face = "bold")
  )
