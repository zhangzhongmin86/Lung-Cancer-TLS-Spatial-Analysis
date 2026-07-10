# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)



#####CD16NK重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "T细胞", "CD16NK", "绘制umap")
df <- fread(file.path(in_dir, "CD16NK_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CCL5_CD16NK" = "#E5D2DD",
  "CD16NK_Cycling"     = "#53A85F",
  "GZMA_CD16NK"     = "#F1BB72",
  "H3-3A_CD16NK"     = "#F3B1A0",
  "HSPA1A_CD16NK"    = "#D6E7A3",
  "JUN_CD16NK"    = "#57C3F3",
  "PRKCH_CD16NK"  = "#476D87",
  "XCL1_CD16NK"   = "#E95C59"
 
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


#####CD56NK重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "T细胞", "CD56NK", "绘制umap")
df <- fread(file.path(in_dir, "CD56NK_obs_umap.csv"))

table(df$celltype_3)



table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CD56NK_Cycling" = "#E5D2DD",
  "GZMH_CD56NK"     = "#53A85F",
  "H3-3A_CD56NK"     = "#F1BB72",
  "IFITM1_CD56NK"     = "#F3B1A0",
  "JUN_CD56NK"    = "#D6E7A3",
  "JUND_CD56NK"    = "#57C3F3",
  "PRKCH_CD56NK"  = "#476D87",
  "XCL1_CD56NK"   = "#E95C59"
  
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

#####NKT重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "T细胞", "NKT", "绘制umap")
df <- fread(file.path(in_dir, "NKT_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "GNB2L1_NKT" = "#E5D2DD",
  "GZMH_NKT"     = "#53A85F",
  "H3-3B_NKT"     = "#F1BB72",
  "JUN_NKT"     = "#F3B1A0",
  "PIP4K2A_NKT"    = "#D6E7A3",
  "PRF1_NKT"    = "#57C3F3",
  "XCL1_NKT"  = "#476D87"
 
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


#####ILT重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "T细胞", "ILT", "绘制umap")
df <- fread(file.path(in_dir, "ILT_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CD74_ILC1-like" = "#E5D2DD",
  "HSP90AA1_ILC3-like"     = "#53A85F",
  "IL32_ILC2-like"     = "#F1BB72",
  "PARP8_ILC1-like"     = "#F3B1A0",
  "TYROBP_ILC3-like"    = "#D6E7A3"
  
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

#####GDT重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "T细胞", "GDT细胞", "绘制umap")
df <- fread(file.path(in_dir, "GDT细胞_obs_umap.csv"))

table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "Anti-tumor_γδ" = "#E5D2DD",
  "Cytotoxic_γδ"     = "#53A85F",
  "gamma_delta_T_Cycling"     = "#F1BB72",
  "Resting_γδ"     = "#F3B1A0",
  "T17_gamma_delta_T"    = "#D6E7A3",
  "Treg_gamma_delta_T"    = "#57C3F3",
  "TRGC2_gamma_delta_T"  = "#476D87",
  "Vδ1_gamma_delta_T"   = "#E95C59",
  "Vδ2_gamma_delta_T"      = "#E59CC4"
  
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











