# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)

######figureS1A######


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

my12colors <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#57C3F3", "#476D87",
  "#E95C59", "#AB3282", "#23452F", "#8C549C", "#585658",
  "#E4C755", "#712820"
)

my36colors <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3",
  "#57C3F3", "#476D87", "#E95C59", "#E59CC4", "#AB3282",
  "#23452F", "#BD956A", "#8C549C", "#585658", "#9FA3A8",
  "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755",
  "#F7F398", "#AA9A59", "#E63863", "#E39A35", "#C1E6F3",
  "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD",
  "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#968175"
)

obs_file <- "/data/beifen/zhongmin/泛癌S/h5ad/189数据/combined_adata_innerobs.csv"

metadata <- read.csv(
  obs_file,
  row.names = 1,
  check.names = TRUE
)

colnames(metadata) <- trimws(colnames(metadata))

required_cols <- c("cancer", "tissue_organ")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop("metadata 中缺少以下列：", paste(missing_cols, collapse = ", "))
}

metadata$cancer <- as.character(metadata$cancer)
metadata$tissue_organ <- as.character(metadata$tissue_organ)

metadata <- metadata %>%
  filter(
    !is.na(cancer),
    cancer != "",
    !is.na(tissue_organ),
    tissue_organ != "",
    tissue_organ != "NA",
    tissue_organ != "nan",
    tissue_organ != "<NA>"
  )

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

cell_counts <- metadata %>%
  group_by(cancer, tissue_organ) %>%
  summarise(cell_number = n(), .groups = "drop")

cell_counts$cancer <- factor(
  cell_counts$cancer,
  levels = cancer_levels_plot
)

tissue_order <- cell_counts %>%
  group_by(tissue_organ) %>%
  summarise(n = sum(cell_number), .groups = "drop") %>%
  arrange(desc(n)) %>%
  pull(tissue_organ)

cell_counts$tissue_organ <- factor(
  cell_counts$tissue_organ,
  levels = tissue_order
)

cell_counts <- cell_counts %>%
  mutate(tissue_numeric = as.numeric(tissue_organ)) %>%
  filter(!is.na(tissue_numeric))

tissues <- levels(cell_counts$tissue_organ)

pal <- if (length(tissues) <= length(my36colors)) {
  my36colors[seq_along(tissues)]
} else {
  rep(my36colors, length.out = length(tissues))
}

tissue_colors <- stats::setNames(pal, tissues)

p_left <- ggplot(cell_counts, aes(x = tissue_numeric, y = cancer)) +
  geom_line(
    aes(group = cancer),
    color = "grey60",
    linewidth = 0.45
  ) +
  geom_point(
    aes(size = cell_number, fill = tissue_organ),
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = tissue_colors,
    drop = FALSE,
    name = "Tissue/Organ"
  ) +
  scale_size(
    range = c(2.8, 11),
    name = "Cell Number"
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(cell_counts$tissue_organ)),
    labels = levels(cell_counts$tissue_organ),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme_classic(
    base_size = base_font_size,
    base_family = base_font_family
  ) +
  theme(
    text = element_text(
      family = base_font_family,
      color = "black"
    ),
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
    axis.line = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.35
    ),
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
  guides(
    size = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  labs(x = NULL, y = NULL)

totals_by_cancer$total_cells_adj <- totals_by_cancer$total_cells + 1

x_min <- max(
  1,
  min(totals_by_cancer$total_cells_adj, na.rm = TRUE) * 0.8
)

x_max <- max(
  totals_by_cancer$total_cells_adj,
  na.rm = TRUE
) * 1.1

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
    xlim = c(x_min, x_max)
  ) +
  theme_classic(
    base_size = base_font_size,
    base_family = base_font_family
  ) +
  theme(
    text = element_text(
      family = base_font_family,
      color = "black"
    ),
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
    axis.line = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.x = element_line(
      color = "black",
      linewidth = 0.35
    ),
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
  rel_widths = c(3.5, 1),
  align = "h",
  axis = "tb"
)

title_plot <- ggdraw() +
  draw_label(
    "Disease Tissue Distribution Analysis",
    fontface = "bold",
    size = main_title_size,
    fontfamily = base_font_family,
    x = 0.5,
    y = 0.70,
    hjust = 0.5
  ) +
  draw_label(
    "Cell counts by tissue/organ and total cells per disease (log scale)",
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
  filename = "tissue_distribution_with_barplot_AI_editable.pdf",
  plot = final_plot,
  device = grDevices::cairo_pdf,
  width = 16,
  height = 14,
  units = "in"
)

ggsave(
  filename = "tissue_distribution_with_barplot_high_resolution.png",
  plot = final_plot,
  width = 16,
  height = 14,
  units = "in",
  dpi = 600
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = "tissue_distribution_with_barplot_AI_editable.svg",
    plot = final_plot,
    device = svglite::svglite,
    width = 16,
    height = 14,
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


######figures1E#######



library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)

obs_file <- "/data/beifen/zhongmin/泛癌S/h5ad/189数据/combined_adata_innerobs.csv"
metadata <- read.csv(obs_file, row.names = 1, check.names = TRUE)

required_cols <- c("group", "celltype_1")
missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop("metadata 中缺少以下列：", paste(missing_cols, collapse = ", "))
}

metadata <- as.data.frame(metadata)

grp_order <- c(
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

grp_order_use <- intersect(grp_order, unique(as.character(metadata$group)))

metadata <- metadata %>%
  dplyr::mutate(
    group = as.character(group),
    celltype_1 = as.character(celltype_1)
  ) %>%
  dplyr::filter(
    !is.na(group),
    group != "",
    !is.na(celltype_1),
    celltype_1 != "",
    group %in% grp_order_use
  ) %>%
  dplyr::mutate(
    group = factor(group, levels = grp_order_use)
  )

cell_color_map <- c(
  "T cells" = '#E5D2DD',
  "Monocyte_macrophage" = '#53A85F',
  "Epithelial cells" = '#F1BB72',
  "B cells" = '#F3B1A0',
  "NK cells" = '#476D87',
  "Endothelial cells" = '#E95C59',
  "Plasma cells" = '#E59CC4',
  "Mast cells" = '#BD956A',
  "Neutrophils" = '#8C549C',
  "DC" = '#23452F',
  "pDC" = '#9FA3A8',
  "Fibroblast cells" = '#57C3F3',
  "Malignant cells" = '#D6E7A3',
  "Alveolar epithelial" = '#AB3282',
  "Hepatic cells" = '#585658'
)

cell_type_order <- names(cell_color_map)

missing_types <- setdiff(unique(metadata$celltype_1), names(cell_color_map))

if (length(missing_types) > 0) {
  message(
    "以下细胞类型未在 cell_color_map 中定义颜色，将从图中移除：",
    paste(missing_types, collapse = ", ")
  )
}

metadata_plot <- metadata %>%
  dplyr::filter(celltype_1 %in% cell_type_order) %>%
  dplyr::mutate(
    celltype_1 = factor(celltype_1, levels = cell_type_order)
  )

label_threshold <- 0.05

comp_df <- metadata_plot %>%
  dplyr::count(group, celltype_1, name = "n") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    perc = n / sum(n),
    label = dplyr::if_else(
      perc >= label_threshold,
      scales::percent(perc, accuracy = 1),
      ""
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    group = factor(group, levels = grp_order_use),
    celltype_1 = factor(celltype_1, levels = cell_type_order)
  )

FONT_FAMILY <- "Arial"

p <- ggplot(comp_df, aes(x = group, y = perc, fill = celltype_1)) +
  geom_col(
    width = 0.92,
    color = "white",
    linewidth = 0.55
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    family = FONT_FAMILY,
    size = 4.8,
    color = "black",
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = cell_color_map,
    breaks = cell_type_order,
    name = "Cell type",
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = grp_order_use,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Fraction"
  ) +
  theme_classic(base_family = FONT_FAMILY) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 0.3),
    axis.line.y = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.18, "cm"),
    axis.text.x = element_text(
      size = 15,
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 15,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 17,
      colour = "black",
      margin = margin(r = 10)
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(
      size = 17,
      face = "bold",
      colour = "black"
    ),
    legend.text = element_text(
      size = 14,
      colour = "black"
    ),
    legend.key.size = unit(0.55, "cm"),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.55, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(12, 18, 12, 12)
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(color = NA)
    )
  )

print(p)

save_dir <- dirname(obs_file)

ggsave(
  filename = file.path(save_dir, "stacked_barplot_topjournal_no_border.pdf"),
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(save_dir, "stacked_barplot_topjournal_no_border.svg"),
  plot = p,
  width = 12.8,
  height = 8.6,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(save_dir, "stacked_barplot_topjournal_no_border.png"),
  plot = p,
  width = 12.8,
  height = 8.6,
  units = "in",
  dpi = 600,
  bg = "white"
)



####figureS1F#####

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

library(ggplot2)
library(dplyr)
library(Seurat)
library(scales)

pbmc = readRDS("/data/beifen/zhongmin/泛癌S/h5ad/189数据/combined_adata_inner_提取10分之一.rds")


celltype_cols <- c(
  "T cells"             = "#E5D2DD",
  "Monocyte_macrophage" = "#53A85F",
  "Epithelial cells"    = "#F1BB72",
  "B cells"             = "#F3B1A0",
  "Fibroblast cells"    = "#57C3F3",
  "NK cells"            = "#476D87",
  "Endothelial cells"   = "#E95C59",
  "Plasma cells"        = "#E59CC4",
  "Alveolar epithelial" = "#AB3282",
  "DC"                  = "#23452F",
  "Mast cells"          = "#BD956A",
  "Neutrophils"         = "#8C549C",
  "Hepatic cells"       = "#585658",
  "pDC"                 = "#9FA3A8",
  "Malignant cells"     = "#D6E7A3"
)

ctype_order <- c(
  "T cells", "NK cells", "B cells", "Plasma cells", "Monocyte_macrophage",
  "DC", "pDC", "Mast cells", "Neutrophils", "Endothelial cells",
  "Fibroblast cells", "Hepatic cells", "Alveolar epithelial", 
  "Epithelial cells", "Malignant cells" 
)

cluster10Marker <- c(
  "CD3D", "CD3E",
  "GNLY", "NKG7",
  "CD79A", "MS4A1",
  "JCHAIN", "IGHA2", "MZB1",
  "CD14", "CD68", "CD163", 
  "CD1C", "CLEC10A", "HLA-DQA1",
  "IRF4", "IRF7", "PLAC8",
  "CPA3", "TPSAB1",
  "ITGAX", "CSF3R", 
  "PECAM1", "CLDN5",
  "COL1A2", "DCN", 
  "PKHD1", "BICC1", "ALB",
  "SFTPB", "SFTPC", 
  "EPCAM", "KRT19"
)



# 

dp <- DotPlot(
  pbmc,
  features       = cluster10Marker,
  group.by       = "celltype_1",  
  assay          = "RNA",
  cluster.idents = TRUE,
  scale          = TRUE,       
  col.min        = -1,
  col.max        = 2,
  dot.min        = 0
)


data <- dp$data %>%
  mutate(
    id            = factor(id, levels = ctype_order),
    features.plot = factor(features.plot, levels = rev(cluster10Marker))  
  )


top_anno <- data %>%
  distinct(id) %>%
  mutate(y_top = length(levels(data$features.plot)) + 1)  



# 
create_bubble_plot_method2_transposed <- function(data, celltype_cols, cluster10Marker) {
  
  # 确保基因的顺序（X轴，从左到右）
  data$features.plot <- factor(data$features.plot, levels = cluster10Marker)
  
  # 准备左侧颜色条数据
  left_anno <- data %>%
    distinct(id) %>%
    mutate(
      x_left = 0, 
      y_numeric = as.numeric(id), 
      color = celltype_cols[as.character(id)]
    )
  
  Fig2e_transposed <- ggplot(data, aes(x = features.plot, y = id)) +
    
    # 绘制左侧的细胞类型颜色侧边栏
    {
      rect_list <- lapply(1:nrow(left_anno), function(i) {
        annotate(
          "rect",
          xmin = -0.4,   
          xmax = 0.4,    
          ymin = left_anno$y_numeric[i] - 0.45, 
          ymax = left_anno$y_numeric[i] + 0.45,
          fill = left_anno$color[i],
          color = NA
        )
      })
      rect_list
    } +
    
    # 气泡矩阵
    geom_point(
      aes(fill = avg.exp.scaled, size = pct.exp),
      shape = 21,
      color = "black",
      stroke = 0.3
    ) +
    
    # 表达量颜色映射
    scale_fill_gradient2(
      low = "#58A4C3",
      mid = "white",
      high = "#E95C59",
      midpoint = 0,
      limits = c(-1, 2.5),
      oob = squish,
      name = "Average\nExpression",
      breaks = c(-1, 0, 1, 2)
    ) +
    
    # 气泡大小映射
    scale_size(
      range = c(0.5, 8),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      name = "Percent\nExpressed"
    ) +
    
    labs(x = NULL, y = NULL) +
    
    # xlim 从 -0.5 开始，给左侧的颜色条留出空间
    coord_cartesian(
      clip = "off",
      xlim = c(-0.5, length(unique(data$features.plot)) + 0.5) 
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      
      # *** 修改点 2: 外部边框变细 ***
      panel.border = element_rect(colour = "black", linewidth = 0.3), # 原来是 1
      
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
      
      # *** 修改点 1: X轴基因文字变大 ***
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 13, # 原来是 10，调大到 13 (可根据需要继续调整)
        face = "italic", 
        color = "black"
      ),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      
      # Y轴：细胞类型 (加粗)
      axis.text.y = element_text(
        face = "bold",
        hjust = 1,
        size = 10,
        color = "black",
        margin = margin(r = 5)
      ),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.6, "cm"),
      legend.position = "right"
    )
  
  return(Fig2e_transposed)
}


data <- DotPlot(pbmc, features = cluster10Marker, group.by = "celltype_1")$data
data$id <- factor(data$id, levels = rev(ctype_order))
data <- data[!is.na(data$id), ]

plot2_transposed <- create_bubble_plot_method2_transposed(data, celltype_cols, cluster10Marker)
print(plot2_transposed)

ggsave("/data/beifen/zhongmin/泛癌S/h5ad/189数据/绘制umap/bubble_plot_method2_transposed.pdf", 
       plot2_transposed, width = 13, height = 9)

