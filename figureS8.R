# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)



 
#####上皮细胞重新绘制umap####

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(ggpubr)
library(ggpubr)
library(ggrepel)
library(grDevices)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(data.table)
library(ggthemes)
library(RColorBrewer)
library(pheatmap)



library(data.table)
library(dplyr)
library(ggplot2)
library(ggrastr)

in_dir <- file.path(PROJECT_ROOT, "转h5ad", "上皮细胞", "绘制umap")
df <- fread(file.path(in_dir, "上皮细胞_obs_umap.csv"))
table(df$celltype_1)
table(df$celltype_2)
table(df$celltype_3)


df <- df %>%
  mutate(
    celltype_2 = as.character(celltype_2),
    celltype_2 = str_trim(celltype_2),
    celltype_2 = paste0(
      str_to_upper(str_sub(celltype_2, 1, 1)),
      str_sub(celltype_2, 2)
    )
  )

table(df$celltype_2)


# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype2_cols <- c(
  "AT0"                       = "#E5D2DD",
  "AT1"                       = "#53A85F",
  "AT2"                       = "#F1BB72",
  "AT2 proliferating"         = "#E95C59",
  "Basal resting"             = "#476D87",
  "Brush cell"                = "#57C3F3",
  "Club"                      = "#D6E7A3",
  "Epithelial cell"           = "#AB3282",
  "Foveolar cell of stomach"  = "#BD956A",
  "Gland cells"               = "#8C549C",
  "Goblet"                    = "#23452F",
  "Hepatocytes"               = "#E59CC4",
  "Hillock-like"              = "#F3B1A0",
  "Ionocyte"                  = "#58A4C3",
  "Malignant cells"           = "#E63863",
  "Multiciliated"             = "#5F3D69",
  "Neuroendocrine"            = "#E4C755",
  "Suprabasal"                = "#585658",
  "Tuft cells"                = "#91D0BE"
)

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_2), names(celltype2_cols))
extra_cols   <- setdiff(names(celltype2_cols), unique(df$celltype_2))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_2 <- factor(df$celltype_2, levels = names(celltype2_cols))

# 3) 计算标签放置位置（各亚群 UMAP 中位数）
# 3) 计算标签放置位置（各亚群 UMAP 中位数）
label_df <- df %>%
  group_by(celltype_2) %>%
  summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

# 2) 手动设置每个标签的偏移量
#    下面这组参数是给 endothelial cells 这张图准备的，
#    你后续如果想更精细，可以继续微调 dx / dy
offset_df <- tribble(
  ~celltype_2,                  ~dx,  ~dy,
  
  "AT0",                         0,   0,
  "AT1",                         0,   0,
  "AT2",                         0,   0.5,
  "AT2 proliferating",           0,   0,
  "Basal resting",               0,   -0.5,
  "Brush cell",                  0,   -0.5,
  "Club",                        0,   0,
  "Epithelial cell",             0,   0,
  "Foveolar cell of stomach",    0,   1,
  "Gland cells",                 1,   -0.5,
  "Goblet",                      0,   0,
  "Hepatocytes",                 0,   -0.3,
  "Hillock-like",                0,   0,
  "Ionocyte",                    0,   0.2,
  "Malignant cells",             0,   0,
  "Multiciliated",               0,   0,
  "Neuroendocrine",              0,   0,
  "Suprabasal",                  0,   0,
  "Tuft cells",                  0,   0
)

# 3) 合并偏移量，得到最终标签位置
label_df <- label_df %>%
  left_join(offset_df, by = "celltype_2") %>%
  mutate(
    dx = ifelse(is.na(dx), 0, dx),
    dy = ifelse(is.na(dy), 0, dy),
    label_x = umap_1 + dx,
    label_y = umap_2 + dy
  )

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_2)) +
  geom_point_rast(
    size = 0.7,
    stroke = 0,
    shape = 16,
    raster.dpi = 600
  ) +
  
  # 白色圆圈画在最终标签位置
  geom_point(
    data = label_df,
    aes(x = label_x, y = label_y),
    color = "black",
    fill = "white",
    size = 7,
    alpha = 0.55,
    shape = 21,
    stroke = 0.4,
    inherit.aes = FALSE
  ) +
  
  # 文字也画在最终标签位置
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = celltype_2),
    colour = "black",
    size = 4.1,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    values = celltype2_cols,
    breaks = names(celltype2_cols),
    name = "Cell type"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 5)
    )
  ) +
  coord_equal(clip = "off") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = unit(5, "mm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(10, 30, 10, 10)
  )

print(p)


#####组织类型####
table(df$tissue_organ)


tissue_cols <- c(
  "Adrenal"      = "#E5D2DD",
  "Bladder"      = "#53A85F",
  "Bone"         = "#F1BB72",
  "Brain"        = "#F3B1A0",
  "breast"       = "#D6E7A3",
  "Esophagus"    = "#57C3F3",
  "Intestine"    = "#476D87",
  "Kidney"       = "#E95C59",
  "Liver"        = "#E59CC4",
  "Lung"         = "#AB3282",
  "Lymph_node"   = "#23452F",
  "Nose"         = "#BD956A",
  "Oral_cavity"  = "#8C549C",
  "Other"        = "#585658",
  "Ovary"        = "#9FA3A8",
  "Stomach"      = "#E0D4CA",
  "Uterus"       = "#5F3D69"
)

# 2) 检查颜色是否完全覆盖
missing_cols <- setdiff(unique(df$tissue_organ), names(tissue_cols))
extra_cols   <- setdiff(names(tissue_cols), unique(df$tissue_organ))

if (length(missing_cols)) {
  message("还没有为这些组织类型指定颜色：", paste(missing_cols, collapse = ", "))
}

if (length(extra_cols)) {
  message("颜色表里这些组织类型在数据中不存在：", paste(extra_cols, collapse = ", "))
}

# 3) 固定图例顺序
df <- df %>%
  mutate(
    tissue_organ = as.character(tissue_organ),
    tissue_organ = factor(tissue_organ, levels = names(tissue_cols))
  )

# 4) 计算标签位置：每个组织类型 UMAP 中位数
label_df <- df %>%
  filter(!is.na(tissue_organ)) %>%
  group_by(tissue_organ) %>%
  summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

# 5) 绘图
p <- ggplot(df, aes(umap_1, umap_2, colour = tissue_organ)) +
  geom_point_rast(
    size = 0.7,
    stroke = 0,
    shape = 16,
    raster.dpi = 600
  ) +
  geom_point(
    data = label_df,
    aes(umap_1, umap_2),
    color = "black",
    fill = "white",
    size = 7,
    alpha = 0.6,
    shape = 21,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = tissue_organ),
    colour = "black",
    size = 3.9,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = tissue_cols,
    breaks = names(tissue_cols),
    name = "Tissue / organ"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 5)
    )
  ) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = unit(5, "mm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank()
  )

print(p)


####celltype_3####
table(df$celltype_3)


tissue_cols <- c(
  "Adrenal"      = "#E5D2DD",
  "Bladder"      = "#53A85F",
  "Bone"         = "#F1BB72",
  "Brain"        = "#F3B1A0",
  "breast"       = "#D6E7A3",
  "Esophagus"    = "#57C3F3",
  "Intestine"    = "#476D87",
  "Kidney"       = "#E95C59",
  "Liver"        = "#E59CC4",
  "Lung"         = "#AB3282",
  "Lymph_node"   = "#23452F",
  "Nose"         = "#BD956A",
  "Oral_cavity"  = "#8C549C",
  "Other"        = "#585658",
  "Ovary"        = "#9FA3A8",
  "Stomach"      = "#E0D4CA",
  "Uterus"       = "#5F3D69"
)

# 2) 检查颜色是否完全覆盖
missing_cols <- setdiff(unique(df$tissue_organ), names(tissue_cols))
extra_cols   <- setdiff(names(tissue_cols), unique(df$tissue_organ))

if (length(missing_cols)) {
  message("还没有为这些组织类型指定颜色：", paste(missing_cols, collapse = ", "))
}

if (length(extra_cols)) {
  message("颜色表里这些组织类型在数据中不存在：", paste(extra_cols, collapse = ", "))
}

# 3) 固定图例顺序
df <- df %>%
  mutate(
    tissue_organ = as.character(tissue_organ),
    tissue_organ = factor(tissue_organ, levels = names(tissue_cols))
  )

# 4) 计算标签位置：每个组织类型 UMAP 中位数
label_df <- df %>%
  filter(!is.na(tissue_organ)) %>%
  group_by(tissue_organ) %>%
  summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

# 5) 绘图
p <- ggplot(df, aes(umap_1, umap_2, colour = tissue_organ)) +
  geom_point_rast(
    size = 0.7,
    stroke = 0,
    shape = 16,
    raster.dpi = 600
  ) +
  geom_point(
    data = label_df,
    aes(umap_1, umap_2),
    color = "black",
    fill = "white",
    size = 7,
    alpha = 0.6,
    shape = 21,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = tissue_organ),
    colour = "black",
    size = 3.9,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = tissue_cols,
    breaks = names(tissue_cols),
    name = "Tissue / organ"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 5)
    )
  ) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = unit(5, "mm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank()
  )

print(p)
