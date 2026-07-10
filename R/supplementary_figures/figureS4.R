# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)

#####肥大细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "髓系细胞", "肥大细胞", "绘制umap")
df <- fread(file.path(in_dir, "肥大细胞_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "GNB2L1_Mast" = "#E5D2DD",
  "HSPA1A_Mast"     = "#53A85F",
  "MCL1_Mast"     = "#F1BB72",
  "MTND4P33_Mast"     = "#F3B1A0",
  "TNF_Mast"    = "#D6E7A3",
  "VEGFA_Mast"    = "#57C3F3"
 
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


#####浆细胞样树突状细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "髓系细胞", "浆细胞样树突状细胞", "绘制umap")
df <- fread(file.path(in_dir, "浆细胞样树突状细胞_obs_umap.csv"))
table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CD74_pDC" = "#E5D2DD",
  "GNB2L1_pDC"     = "#53A85F",
  "GPR183_pDC"     = "#F1BB72",
  "HSPA1A_pDC"     = "#F3B1A0",
  "ISG15_pDC"    = "#D6E7A3",
  "S100A8_pDC"    = "#57C3F3"
 
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


#####中性粒细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "髓系细胞", "中性粒细胞", "绘制umap")
df <- fread(file.path(in_dir, "中性粒细胞_obs_umap.csv"))
table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CXCL8_CD74_Neutrophils" = "#E5D2DD",
  "CXCL8_Neutrophils"     = "#53A85F",
  "CXCR2_Neutrophils"     = "#F1BB72",
  "FCGR3B_Neutrophils"     = "#F3B1A0",
  "IFIT1_Neutrophils"    = "#D6E7A3",
  "MALAT1_Neutrophils"    = "#57C3F3",
  "MME_Neutrophils"  = "#476D87",
  "MMP9_Neutrophils"   = "#E95C59",
  "NFKBIZ_Neutrophils"      = "#E59CC4"
  
)

# （可选）检查是否有未覆盖或多余的类别
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
  ~celltype_3,                  ~dx,   ~dy,
  
  "CXCL8_CD74_Neutrophils",     -1.4,  0.5,
  "CXCL8_Neutrophils",           1.0,  0.4,
  "CXCR2_Neutrophils",          -0.6,  0.0,
  "FCGR3B_Neutrophils",          0.3,  0.2,
  "IFIT1_Neutrophils",           0.0,  0.9,
  "MALAT1_Neutrophils",         -1.0, -0.2,
  "MME_Neutrophils",             0.0, -0.0,
  "MMP9_Neutrophils",            0.8, -0.1,
  "NFKBIZ_Neutrophils",          0.3, -3.5
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


