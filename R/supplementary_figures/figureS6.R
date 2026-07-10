# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)




#####CD4T重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "转h5ad", "T细胞", "CD4", "绘制umap")
df <- fread(file.path(in_dir, "CD4_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CD4 T_Cycling"             = "#476D87",
  "CD4 TCM"                   = "#E5D2DD",
  "CD4 TEM"                   = "#53A85F",
  "CD4 Trm"                   = "#F1BB72",
  "COX2_CD4 TEM"              = "#E95C59",
  "ISG15_Treg"                = "#AB3282",
  "MT1X_Naive CD4 T"          = "#57C3F3",
  "Naive CD4 T"               = "#D6E7A3",
  "Pre_Th"                    = "#8C549C",
  "Stress response CD4 TCM"   = "#E63863",
  "Stress response CD4 TEM"   = "#E39A35",
  "Tfh"                       = "#58A4C3",
  "Th17"                      = "#23452F",
  "TMSB4X_CD4 T"              = "#BD956A",
  "Treg"                      = "#585658"
  
)

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_3), names(celltype3_cols))
extra_cols   <- setdiff(names(celltype3_cols), unique(df$celltype_3))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_3 <- factor(df$celltype_3, levels = names(celltype3_cols))

# 3) 计算标签放置位置（各亚群 UMAP 中位数）
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
  
  "CD4 T_Cycling",               0.0,   0.9,
  "CD4 TCM",                    -0.7,  -1.9,
  "CD4 TEM",                    -1.2,  -0.1,
  "CD4 Trm",                    -0.8,   0.2,
  "COX2_CD4 TEM",                0.8,   0.2,
  "ISG15_Treg",                  0.6,   0.9,
  "MT1X_Naive CD4 T",            0.3,   0.0,
  "Naive CD4 T",                -0.4,  -2.3,
  "Pre_Th",                      1.5,   0.6,
  "Stress response CD4 TCM",    -0.1,   0.0,
  "Stress response CD4 TEM",    -0.7,   0.4,
  "Tfh",                         0.0,   1.2,
  "Th17",                        0.0,   0.0,
  "TMSB4X_CD4 T",                1.1,  -0.4,
  "Treg",                        0.4,   0.8
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











#####CD8T重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "转h5ad", "T细胞", "CD8", "绘制umap")
df <- fread(file.path(in_dir, "CD8_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CD8 TCM"                            = "#E5D2DD",
  "CD8 TEM"                            = "#53A85F",
  "CD8 Temra"                          = "#F1BB72",
  "CD8 Trm"                            = "#F3B1A0",
  "Exhausted CD8 TEM"                  = "#E95C59",
  "Exhausted CD8_Cycling"              = "#476D87",
  "Exhausted_Stress response CD8 TEM"  = "#AB3282",
  "ISG15 CD8 TEFF"                     = "#57C3F3",
  "MAIT"                               = "#23452F",
  "Naive CD8 T"                        = "#D6E7A3",
  "Naive CD8_Cycling"                  = "#E63863",
  "Precursor exhausted CD8 TCM"        = "#8C549C",
  "Precursor exhausted CD8 TEM"        = "#BD956A",
  "STAT4_CD8"                          = "#58A4C3",
  "TMSB4X_CD8 T"                       = "#585658"
  
)

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_3), names(celltype3_cols))
extra_cols   <- setdiff(names(celltype3_cols), unique(df$celltype_3))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_3 <- factor(df$celltype_3, levels = names(celltype3_cols))

# 3) 计算标签放置位置（各亚群 UMAP 中位数）
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
  ~celltype_3,                         ~dx,   ~dy,
  
  "CD8 TCM",                           -2.5,   0.2,
  "CD8 TEM",                           -0.6,   -0.7,
  "CD8 Temra",                         -1.1,  -0.2,
  "CD8 Trm",                            1.0,   0.3,
  "Exhausted CD8 TEM",                 -0.8,  -0.3,
  "Exhausted CD8_Cycling",              0.4,  -0.9,
  "Exhausted_Stress response CD8 TEM",  0.5,  -1.0,
  "ISG15 CD8 TEFF",                    -0.9,  -3.7,
  "MAIT",                               1.8,   1.2,
  "Naive CD8 T",                       -0.2,   1.7,
  "Naive CD8_Cycling",                 -1.1,   0.5,
  "Precursor exhausted CD8 TCM",        1.6,  -4.1,
  "Precursor exhausted CD8 TEM",        0.6,  -2.3,
  "STAT4_CD8",                         -1.0,   0.6,
  "TMSB4X_CD8 T",                       0.1,   0.9
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









