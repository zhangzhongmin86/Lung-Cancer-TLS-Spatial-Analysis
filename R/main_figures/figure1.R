# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)


library(Seurat)
library(Matrix)
library(readr)
library("moduleColor")
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
#设置图片输出目录


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)



######figure1C######
library(dplyr)
library(ggplot2)
library(cowplot)

base_font_family <- "Arial"
base_font_size <- 14

axis_x_size <- 13
axis_y_size <- 13
legend_text_size <- 12
legend_title_size <- 13
right_axis_x_size <- 12
right_axis_title_size <- 13
main_title_size <- 22
subtitle_size <- 15

theme_set(
  theme_classic(base_size = base_font_size, base_family = base_font_family) +
    theme(
      text = element_text(family = base_font_family, color = "black"),
      axis.text = element_text(family = base_font_family, color = "black"),
      axis.title = element_text(family = base_font_family, color = "black"),
      legend.text = element_text(
        family = base_font_family,
        size = legend_text_size,
        color = "black"
      ),
      legend.title = element_text(
        family = base_font_family,
        size = legend_title_size,
        face = "bold",
        color = "black"
      )
    )
)

obs_file <- "/data/beifen/zhongmin/泛癌S/h5ad/189数据/combined_adata_innerobs.csv"
metadata <- read.csv(obs_file, row.names = 1, check.names = TRUE)
colnames(metadata) <- trimws(colnames(metadata))

required_cols <- c("cancer", "group")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop("metadata 中缺少以下列：", paste(missing_cols, collapse = ", "))
}

metadata$cancer <- as.character(metadata$cancer)
metadata$group <- as.character(metadata$group)

totals_by_cancer <- metadata %>%
  group_by(cancer) %>%
  summarise(total_cells = n(), .groups = "drop") %>%
  arrange(total_cells)

cancer_levels_asc <- totals_by_cancer$cancer
cancer_levels_plot <- rev(cancer_levels_asc)

totals_by_cancer$cancer <- factor(
  totals_by_cancer$cancer,
  levels = cancer_levels_plot
)

cell_counts_group <- metadata %>%
  group_by(cancer, group) %>%
  summarise(cell_number = n(), .groups = "drop")

group_levels <- c(
  "Autoimmune diseases",
  "Inflammation",
  "Normal",
  "Normal adjacent tissue",
  "Precancerous condition",
  "Tumor",
  "Tumor_metastasis",
  "Tumor_blood",
  "Blood_Tumor"
)

unexpected_groups <- setdiff(unique(cell_counts_group$group), group_levels)

if (length(unexpected_groups) > 0) {
  message("以下 group 不在预设 group_levels 中，绘图时会被设为 NA：")
  message(paste(unexpected_groups, collapse = ", "))
}

cell_counts_group$cancer <- factor(
  cell_counts_group$cancer,
  levels = cancer_levels_plot
)

cell_counts_group$group <- factor(
  cell_counts_group$group,
  levels = group_levels
)

cell_counts_group <- cell_counts_group %>%
  mutate(group_numeric = as.numeric(group)) %>%
  filter(!is.na(group_numeric))

group_colors <- c(
  "Autoimmune diseases" = "#F1BB72",
  "Inflammation" = "#D6E7A3",
  "Normal" = "#E5D2DD",
  "Normal adjacent tissue" = "#57C3F3",
  "Precancerous condition" = "#E59CC4",
  "Tumor" = "#53A85F",
  "Tumor_metastasis" = "#F3B1A0",
  "Tumor_blood" = "#476D87",
  "Blood_Tumor" = "#E95C59"
)

p_left <- ggplot(cell_counts_group, aes(x = group_numeric, y = cancer)) +
  geom_line(
    aes(group = cancer),
    color = "grey60",
    linewidth = 0.45
  ) +
  geom_point(
    aes(size = cell_number, fill = group),
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE,
    name = "Group"
  ) +
  scale_size(
    range = c(2.8, 11),
    name = "Cell Number"
  ) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme_classic(
    base_size = base_font_size,
    base_family = base_font_family
  ) +
  theme(
    text = element_text(family = base_font_family, color = "black"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = axis_x_size,
      family = base_font_family,
      color = "black"
    ),
    axis.text.y = element_text(
      size = axis_y_size,
      face = "bold",
      family = base_font_family,
      color = "black"
    ),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    legend.position = "left",
    legend.text = element_text(
      size = legend_text_size,
      family = base_font_family,
      color = "black"
    ),
    legend.title = element_text(
      size = legend_title_size,
      face = "bold",
      family = base_font_family,
      color = "black"
    ),
    plot.margin = margin(10, 5, 10, 10)
  ) +
  labs(x = NULL, y = NULL)

totals_by_cancer$total_cells_adj <- totals_by_cancer$total_cells + 1

p_right <- ggplot(totals_by_cancer, aes(y = cancer, x = total_cells_adj)) +
  geom_col(
    fill = "steelblue",
    width = 0.7,
    alpha = 0.8
  ) +
  scale_x_log10(
    breaks = c(1e3, 1e4, 1e5, 1e6),
    labels = c("1K", "10K", "100K", "1M"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  coord_cartesian(
    xlim = c(800, max(totals_by_cancer$total_cells_adj, na.rm = TRUE) * 1.1)
  ) +
  theme_classic(
    base_size = base_font_size,
    base_family = base_font_family
  ) +
  theme(
    text = element_text(family = base_font_family, color = "black"),
    axis.title.x = element_text(
      size = right_axis_title_size,
      family = base_font_family,
      color = "black"
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      size = right_axis_x_size,
      angle = 45,
      hjust = 1,
      vjust = 1,
      family = base_font_family,
      color = "black"
    ),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linewidth = 0.35),
    plot.margin = margin(10, 10, 10, 0),
    panel.grid.major.x = element_line(
      color = "grey90",
      linetype = "dashed",
      linewidth = 0.35
    ),
    panel.grid.minor.x = element_line(
      color = "grey95",
      linetype = "dotted",
      linewidth = 0.25
    )
  ) +
  labs(x = "Total Cells (log scale)", y = NULL)

combined_plot <- plot_grid(
  p_left,
  p_right,
  nrow = 1,
  rel_widths = c(3, 1),
  align = "h",
  axis = "tb"
)

title_plot <- ggdraw() +
  draw_label(
    "Disease Cell Distribution Analysis",
    fontface = "bold",
    size = main_title_size,
    fontfamily = base_font_family,
    x = 0.5,
    y = 0.70,
    hjust = 0.5
  ) +
  draw_label(
    "Cell counts by group and total cells per disease (log scale)",
    size = subtitle_size,
    fontfamily = base_font_family,
    x = 0.5,
    y = 0.25,
    hjust = 0.5
  )

final_plot <- plot_grid(
  title_plot,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.10, 1)
)

print(final_plot)

ggsave(
  filename = "group_distribution_with_barplot_AI_editable.pdf",
  plot = final_plot,
  device = cairo_pdf,
  width = 10,
  height = 13,
  units = "in"
)

ggsave(
  filename = "group_distribution_with_barplot_high_resolution.png",
  plot = final_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 600
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = "group_distribution_with_barplot_AI_editable.svg",
    plot = final_plot,
    device = svglite::svglite,
    width = 14,
    height = 12,
    units = "in"
  )
} else {
  message("未安装 svglite 包，跳过 SVG 保存。PDF 已保存，可直接用于 Adobe Illustrator。")
}

cat("\n细胞数量统计:\n")
cat("最小值:", min(totals_by_cancer$total_cells, na.rm = TRUE), "\n")
cat("最大值:", max(totals_by_cancer$total_cells, na.rm = TRUE), "\n")
cat("中位数:", median(totals_by_cancer$total_cells, na.rm = TRUE), "\n")
cat(
  "数据跨越的数量级:",
  round(
    log10(
      max(totals_by_cancer$total_cells, na.rm = TRUE) /
        min(totals_by_cancer$total_cells, na.rm = TRUE)
    ),
    1
  ),
  "\n\n"
)



