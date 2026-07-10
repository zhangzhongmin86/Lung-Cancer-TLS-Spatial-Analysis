# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)



#####成纤维细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "成纤维细胞", "绘制umap")
df <- fread(file.path(in_dir, "成纤维细胞_obs_umap.csv"))

# 1) 清洗类别名，修正 “MINOS1 _Fib” 这类带空格的情况
df <- df %>%
  mutate(
    celltype_3 = trimws(as.character(celltype_3)),
    celltype_3 = sub("^MINOS1\\s+_Fib$", "MINOS1_Fib", celltype_3)  # 如需可再加别名修正
  )

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "ADAMDEC1_Fib" = "#E5D2DD",
  "CCL2_Fib"     = "#53A85F",
  "CD74_Fib"     = "#F1BB72",
  "DLX5_Fib"     = "#F3B1A0",
  "DOCK4_Fib"    = "#D6E7A3",
  "FBXL7_Fib"    = "#57C3F3",
  "Fib_Cycling"  = "#476D87",
  "GNB2L1_Fib"   = "#E95C59",
  "GSN_Fib"      = "#E59CC4",
  "MINOS1_Fib"   = "#AB3282",
  "MMP1_Fib"     = "#23452F",
  "MMP11_Fib"    = "#BD956A",
  "MYH11_SMC"    = "#8C549C",
  "Pericyte"     = "#585658",
  "PI16_Fib"     = "#9FA3A8",
  "RGCC_Fib"     = "#E0D4CA",
  "SOX6_Fib"     = "#5F3D69"
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
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_3)) +
  geom_point_rast(size = 0.7, stroke = 0, shape = 16, raster.dpi = 600) +
  geom_point(
    data = label_df, aes(umap_1, umap_2),
    color = "black", fill = "white", size = 7, alpha = 0.6, shape = 21, inherit.aes = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = celltype_3),
    colour = "black",
    size = 3.9,
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

ggsave("~/daima/文章整理代码/补充图S2/Figure1A_fromCSV_celltype3_fixed_colors.pdf",
       p, width = 9, height = 6)


table(df$celltype_2)
# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype2_cols <- c(
  "Fib_alveolar" = "#E5D2DD",
  "Fib_Cycling"     = "#53A85F",
  "Fib_epithelial crypt"     = "#F1BB72",
  "Fib_i"     = "#F3B1A0",
  "Fib_lamina propria"    = "#D6E7A3",
  "Fib_my"    = "#57C3F3",
  "Fib_progenitor-like"  = "#476D87",
  "Pericyte"   = "#E95C59",
  "SMC"      = "#E59CC4"
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



#####内皮细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "内皮细胞", "绘制umap")
df <- fread(file.path(in_dir, "内皮细胞_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "Arterial ECs" = "#E5D2DD",
  "ECs_Cycling"     = "#53A85F",
  "EDNRB_Capillary ECs"     = "#F1BB72",
  "EPAS1_Capillary ECs"     = "#F3B1A0",
  "H3-3B_Tip Cells"    = "#D6E7A3",
  "IGFBP5_Capillary ECs"    = "#57C3F3",
  "LDB2_Capillary ECs"  = "#476D87",
  "Lymphatic ECs"   = "#E95C59",
  "MT2A_Capillary ECs"      = "#E59CC4",
  "PLVAP_Capillary ECs"   = "#AB3282",
  "SEPW1_Capillary ECs"     = "#23452F",
  "SPARCL1_Capillary ECs"    = "#BD956A",
  "Tip Cells"    = "#8C549C",
  "Venous ECs"     = "#585658"
  
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
  ~celltype_3,                ~dx,   ~dy,
  "Arterial ECs",              1.3, -0.1,
  "ECs_Cycling",               0.1,  0.7,
  "EDNRB_Capillary ECs",      -1.1,  0.3,
  "EPAS1_Capillary ECs",       0.8,  0.6,
  "H3-3B_Tip Cells",          -0.8,  0.2,
  "IGFBP5_Capillary ECs",      1.0,  0.1,
  "LDB2_Capillary ECs",       -3.0, -0.7,
  "Lymphatic ECs",            -1.0,  0.3,
  "MT2A_Capillary ECs",        0.9, -1.2,
  "PLVAP_Capillary ECs",      -0.2, -0.2,
  "SEPW1_Capillary ECs",      -1.0, -0.2,
  "SPARCL1_Capillary ECs",    -1.0, -0.5,
  "Tip Cells",                 1.7,  0.1,
  "Venous ECs",                0.2, -0.9
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
  "Arterial ECs" = "#E5D2DD",
  "Capillary ECs"     = "#53A85F",
  "ECs_Cycling"     = "#F1BB72",
  "Lymphatic ECs"     = "#F3B1A0",
  "Tip Cells"    = "#D6E7A3",
  "Venous ECs"    = "#57C3F3"
  
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










