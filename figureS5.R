# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)




#####B细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "B细胞", "绘制umap")
df <- fread(file.path(in_dir, "B细胞_obs_umap.csv"))
table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "CXCR4_B memory" = "#E5D2DD",
  "FCHSD2_B memory"     = "#53A85F",
  "GNB2L1_B memory"     = "#F1BB72",
  "H3-3A_B memory"     = "#F3B1A0",
  "HSPA1A_B memory"    = "#D6E7A3",
  "IFIT3_B memory"    = "#57C3F3",
  "S100A9_B memory"  = "#476D87",
  "SOX4_B_cell"   = "#E95C59"
 
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


#####浆细胞重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "浆细胞", "绘制umap")
df <- fread(file.path(in_dir, "浆细胞_obs_umap.csv"))
table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "HSPA1A_plasma cells" = "#E5D2DD",
  "IGHA1_plasma cells"     = "#53A85F",
  "IGHG3_plasma cells"     = "#F1BB72",
  "IGJ_plasma cells"     = "#F3B1A0",
  "IGKV1-12_plasma cells"    = "#D6E7A3",
  "IGLC3_plasma cells"    = "#57C3F3",
  "IGLL5_plasma cells"  = "#476D87",
  "IGLV2-14_plasma cells"   = "#E95C59"
  
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



#####浆母重新绘制umap####

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

in_dir <- file.path(PROJECT_ROOT, "泛癌S", "h5ad", "浆母细胞", "绘制umap")
df <- fread(file.path(in_dir, "浆母细胞_obs_umap.csv"))
table(df$celltype_3)

# 2) 为 celltype_3 明确指定“名字=颜色”的映射（你可以改成自己喜欢的颜色）
celltype3_cols <- c(
  "GCB" = "#E5D2DD",
  "Plasmablast"     = "#53A85F"
  
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






