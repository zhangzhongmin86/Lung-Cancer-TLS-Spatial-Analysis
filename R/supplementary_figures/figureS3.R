# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)



#####单核巨噬细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "整理", "髓系细胞", "单核巨噬细胞", "绘制umap")
df <- fread(file.path(in_dir, "单核巨噬细胞_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  
    "APOC1_Macrophage"        = "#E5D2DD",
    "C1QC_Macrophage"         = "#53A85F",
    "CCSER1_Macrophage"       = "#F1BB72",
    "CD14CD16 monocyte"       = "#F3B1A0",
    "CTSL_Macrophage"         = "#D6E7A3",
    "DNASE1L3_Macrophage"     = "#57C3F3",
    "GABARAP_CD14 monocyte"   = "#476D87",
    "GNB2L1_Macrophage"       = "#E95C59",
    "H3-3A_CD14 monocyte"     = "#E59CC4",
    "HSPA1A_Macrophage"       = "#AB3282",
    "IFI30_CD14 monocyte"     = '#B53E2B',
    "IL1B_CD14 monocyte"      = "#BD956A",
    "ISG15_Macrophage"        =  '#3A6963',
    "LRMDA_Macrophage"        = "#585658",
    "LST1_CD16 monocyte"      = "#9FA3A8",
    "Macrophage_Cycling"      = "#E0D4CA",
    "MT2A_Macrophage"         = "#5F3D69",
    "PDZK1IP1_Macrophage"     = '#E39A35',
    "S100A8_CD14 monocyte"    = "#58A4C3",
    "SPARCL1_Macrophage"      = "#E4C755"
  )
  


# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_3), names(celltype3_cols))
extra_cols   <- setdiff(names(celltype3_cols), unique(df$celltype_3))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_3 <- factor(df$celltype_3, levels = names(celltype3_cols))


# 3) 计算标签放置位置（各亚群 UMAP 中位数）
label_df <- df %>%
  group_by(celltype_3) %>%
  summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

# 2) 手动设置每个标签的偏移量
#    下面这组参数是给 endothelial cells 这张图准备的，
#    你后续如果想更精细，可以继续微调 dx / dy
offset_df <- tribble(
  ~celltype_3,                 ~dx,   ~dy,
  
  "APOC1_Macrophage",          -1.0,  0.3,
  "C1QC_Macrophage",           -0.7,  0.5,
  "CCSER1_Macrophage",          1.2,  0.6,
  "CD14CD16 monocyte",         -0.1, -0.2,
  "CTSL_Macrophage",            0.5, -0.5,
  "DNASE1L3_Macrophage",       -1.0,  0.8,
  "GABARAP_CD14 monocyte",      2.4,  0.2,
  "GNB2L1_Macrophage",        -1.5,  4.0,
  "H3-3A_CD14 monocyte",       -0.4, -0.2,
  "HSPA1A_Macrophage",         -3.8,  2.8,
  "IFI30_CD14 monocyte",        0.3,  0.5,
  "IL1B_CD14 monocyte",         0.9, -0.1,
  "ISG15_Macrophage",           0.6,  0.3,
  "LRMDA_Macrophage",           0.1,  0.3,
  "LST1_CD16 monocyte",        -1.5, -0.4,
  "Macrophage_Cycling",        -0.1,  1.0,
  "MT2A_Macrophage",           -0.1,  0.9,
  "PDZK1IP1_Macrophage",       -0.1,  0.3,
  "S100A8_CD14 monocyte",       0.5, -0.5,
  "SPARCL1_Macrophage",        -0.1,  4.2
)

# 3) 合并偏移量，得到最终标签位置
label_df <- label_df %>%
  left_join(offset_df, by = "celltype_3") %>%
  mutate(
    dx = ifelse(is.na(dx), 0, dx),
    dy = ifelse(is.na(dy), 0, dy),
    label_x = umap_1 + dx,
    label_y = umap_2 + dy
  )

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_3)) +
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
    aes(x = label_x, y = label_y, label = celltype_3),
    colour = "black",
    size = 4.1,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    values = celltype3_cols,
    breaks = names(celltype3_cols),
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

ggsave(file.path(in_dir, "Figure1A_fromCSV_celltype3_fixed_colors.pdf"),
       p, width = 9, height = 6)


table(df$celltype_2)
# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype2_cols <- c(
  "Alveolar macrophage" = "#E5D2DD",
  "CD14 monocyte"     = "#53A85F",
  "CD14CD16 monocyte"     = "#F1BB72",
  "CD16 monocyte"     = "#F3B1A0",
  
  "Macrophage"    = "#57C3F3"
  
)

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_2), names(celltype2_cols))
extra_cols   <- setdiff(names(celltype2_cols), unique(df$celltype_2))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_2 <- factor(df$celltype_2, levels = names(celltype2_cols))

# 3) 计算标签放置位置（各亚群 UMAP 中位数）
label_df <- df %>%
  group_by(celltype_2) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_2)) +
  geom_point_rast(size = 0.7, stroke = 0, shape = 16, raster.dpi = 600) +
  geom_point(
    data = label_df, aes(umap_1, umap_2),
    color = "black", fill = "white", size = 7, alpha = 0.6, shape = 21, inherit.aes = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = celltype_2),
    colour = "black",
    size = 3.9,
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



#####树突状细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "髓系细胞", "树突状细胞", "绘制umap")
df <- fread(file.path(in_dir, "树突状细胞_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "DC1" = "#E5D2DD",
  "LAMP3_DC"     = "#53A85F",
  "Langerhans"     = "#F1BB72",
  "ATP5F1E_DC2"     = "#F3B1A0",
  "CD14_DC2"    = "#D6E7A3",
  "HSPA1A_DC2"    = "#57C3F3",
  "IL1B_DC2"  = "#476D87",
  "IL6ST_DC2"   = "#E95C59",
  "ISG15_DC2"      = "#E59CC4",
  "DC_Cycling"   = "#AB3282"
 
)

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_3), names(celltype3_cols))
extra_cols   <- setdiff(names(celltype3_cols), unique(df$celltype_3))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

df$celltype_3 <- factor(df$celltype_3, levels = names(celltype3_cols))


# 3) 计算标签放置位置（各亚群 UMAP 中位数）
label_df <- df %>%
  group_by(celltype_3) %>%
  summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

# 2) 手动设置每个标签的偏移量
#    下面这组参数是给 endothelial cells 这张图准备的，
#    你后续如果想更精细，可以继续微调 dx / dy
offset_df <- tribble(
  ~celltype_3,        ~dx,   ~dy,
  
  "DC1",              0.2,  -1.1,
  "LAMP3_DC",         0.1,   0.8,
  "Langerhans",       1.0,  -0.1,
  "ATP5F1E_DC2",      0.0,  -0.8,
  "CD14_DC2",        -0.8,  -0.3,
  "HSPA1A_DC2",       0.6,   0.3,
  "IL1B_DC2",        -0.4,   0.4,
  "IL6ST_DC2",       -0.8,   0.5,
  "ISG15_DC2",        1.4,   1.3,
  "DC_Cycling",       0.9,   0.5
)

# 3) 合并偏移量，得到最终标签位置
label_df <- label_df %>%
  left_join(offset_df, by = "celltype_3") %>%
  mutate(
    dx = ifelse(is.na(dx), 0, dx),
    dy = ifelse(is.na(dy), 0, dy),
    label_x = umap_1 + dx,
    label_y = umap_2 + dy
  )

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_3)) +
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
    aes(x = label_x, y = label_y, label = celltype_3),
    colour = "black",
    size = 4.1,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    values = celltype3_cols,
    breaks = names(celltype3_cols),
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
ggsave(file.path(in_dir, "Figure1A_fromCSV_celltype3_fixed_colors.pdf"),
       p, width = 9, height = 6)


table(df$celltype_2)
# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype2_cols <- c(
  "DC1" = "#E5D2DD",
  "LAMP3_DC"     = "#53A85F",
  "Langerhans"     = "#F1BB72",
  "DC2"     = "#F3B1A0",
 
  "DC_Cycling"   = "#AB3282"
)

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_2), names(celltype2_cols))
extra_cols   <- setdiff(names(celltype2_cols), unique(df$celltype_2))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_2 <- factor(df$celltype_2, levels = names(celltype2_cols))

# 3) 计算标签放置位置（各亚群 UMAP 中位数）
label_df <- df %>%
  group_by(celltype_2) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_2)) +
  geom_point_rast(size = 0.7, stroke = 0, shape = 16, raster.dpi = 600) +
  geom_point(
    data = label_df, aes(umap_1, umap_2),
    color = "black", fill = "white", size = 7, alpha = 0.6, shape = 21, inherit.aes = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = celltype_2),
    colour = "black",
    size = 3.9,
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










