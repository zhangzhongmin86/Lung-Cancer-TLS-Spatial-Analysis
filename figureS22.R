# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)

####绘制umap#####


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

in_dir <- "/data/beifen/zhongmin/slide-tag/H5AD格式数据/导出umap"
df <- fread(file.path(in_dir, "slide_tag_lung_obs_umap.csv"))
table(df$celltype_3)
table(df$celltype_3_ZZM)



df$celltype_3_ZZM <- as.character(df$celltype_3_ZZM)

# 把 SMC 和 Pericyte 改成 Fibroblasts
# df$celltype_3_ZZM[df$celltype_3_ZZM %in% c("SMC", "Pericyte")] <- "Fibroblasts"
table(df$celltype_3_ZZM)
colnames(df)
# 你的 36 色库（按出现顺序）
my36colors <- c(
  '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
  '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
  '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
  '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
  '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
  '#968175'
)


celltype3_cols <- c(
  "Arterial ECs"   = "#E5D2DD",
  "AT1"            = "#53A85F",
  "AT2"            = "#F1BB72",
  "B cells"        = "#F3B1A0",
  "Capillary ECs"  = "#D6E7A3",
  "CD4 T cells"    = "#57C3F3",
  "CD8 T cells"    = "#476D87",
  "Malignant cells"= "#E95C59",
  "DC2"            = "#E59CC4",
  "Fibroblasts"    = "#AB3282",
  "LAMP3_DC"       = "#23452F",
  "Lymphatic ECs"  = "#BD956A",
  "Macrophage"     = "#8C549C",
  # "M2-like"     = '#968175',
  
  "DC1"            = "#585658",
  "Mast cells"     = "#9FA3A8",
  "Multiciliated"  = "#E0D4CA",
  "NK cells"       = "#5F3D69",
  
  "Plasma cells"   = "#58A4C3",
  
  "Tip cells"      = '#B53E2B',#"#F7F398",
  "Treg"           = "#AA9A59",
  "Venous ECs"     = '#E39A35',
  "Pericyte"       = "#C5DEBA",
  "SMC"            = "#E4C755",
  "CD4 Trm"        ='#968175',
  "Tfh"            ='#712820',
  "CD4 T_Cycling" = '#625D9E',
  "CD8 Trm"      =  '#F7F398',
  "B naive"      = '#625D9E',
  "GCB"           = '#68A180',
  "FDC"          =   '#3A6963'
  
)

# 如果后续还有新细胞类型，继续往下 my36colors[27] ... 即可

# （可选）检查是否有未覆盖或多余的类别
missing_cols <- setdiff(unique(df$celltype_3_ZZM), names(celltype3_cols))
extra_cols   <- setdiff(names(celltype3_cols), unique(df$celltype_3_ZZM))
if (length(missing_cols)) message("还没有为这些类别指定颜色： ", paste(missing_cols, collapse=", "))
if (length(extra_cols))   message("颜色表里这些类别在数据中不存在： ", paste(extra_cols, collapse=", "))

# 固定图例顺序（用你定义的顺序）
df$celltype_3_ZZM <- factor(df$celltype_3_ZZM, levels = names(celltype3_cols))

# 3) 计算标签放置位置（各亚群 UMAP 中位数）
label_df <- df %>%
  group_by(celltype_3_ZZM) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

# 4) 画图
p <- ggplot(df, aes(umap_1, umap_2, colour = celltype_3_ZZM)) +
  geom_point_rast(size = 0.5, stroke = 0, shape = 16, raster.dpi = 600) +
  geom_point(
    data = label_df, aes(umap_1, umap_2),
    color = "black", fill = "white", size = 4.5, alpha = 0.6, shape = 21, inherit.aes = FALSE
  ) +
  
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = celltype_3_ZZM),
    colour = "black",
    size = 2.5,
    inherit.aes = FALSE
  )+
  
  scale_colour_manual(
    values = celltype3_cols,
    breaks = names(celltype3_cols),
    name   = "Cell type (fibro subclusters)"
  ) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
  )

print(p)

# ggsave(file.path(in_dir, "Figure1A_fromCSV_celltype3_fixed_colors.pdf"),
#        p, width = 12, height = 9)


table(df$sample)

#按照样本着色

# ===== 按“样本 sample”着色，无样本名标注；使用你给的 my36colors =====
library(ggplot2)
library(dplyr)



table(df$sample)

my36colors <- c(
  '#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87',
  '#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658',
  '#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755','#F7F398',
  '#AA9A59','#E63863','#E39A35','#C1E6F3','#6778AE','#91D0BE','#B53E2B',
  '#712820','#DCC1DD','#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963',
  '#968175'
)

df <- df[
  !is.na(df$sample) &
    !is.na(df$umap_1) &
    !is.na(df$umap_2),
]

samps <- sort(unique(df$sample))

df$sample <- factor(df$sample, levels = samps)

sample_cols <- setNames(
  rep_len(my36colors, length(samps)),
  samps
)

p <- ggplot(df, aes(x = umap_1, y = umap_2, colour = sample)) +
  ggrastr::geom_point_rast(
    size = 0.6,
    stroke = 0,
    shape = 16,
    alpha = 0.9,
    raster.dpi = 300
  ) +
  scale_colour_manual(
    values = sample_cols,
    breaks = names(sample_cols),
    name = "Sample"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 5, alpha = 1)
    )
  ) +
  coord_equal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = grid::unit(5, "mm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    
    panel.border = element_blank(),
    axis.line = element_blank()
  )

print(p)



#####按亚型着色#####

library(ggplot2)
library(dplyr)


table(df$subtype)

subtype_cols <- c(
  "Normal adjacent tissue" = "#57C3F3",
  "AIS"                    = "#F1BB72",
  "MIA"                    = "#E59CC4",
  "IA"                     = "#53A85F",
  "COPD"                   = "#E95C59",
  "Granulomatosis"         = "#476D87"
)

missing_cols <- setdiff(unique(df$subtype), names(subtype_cols))
extra_cols <- setdiff(names(subtype_cols), unique(df$subtype))

if (length(missing_cols)) {
  message("还没有为这些 subtype 指定颜色：", paste(missing_cols, collapse = ", "))
}
if (length(extra_cols)) {
  message("颜色表里这些 subtype 在数据中不存在：", paste(extra_cols, collapse = ", "))
}

df <- df[
  !is.na(df$subtype) &
    !is.na(df$umap_1) &
    !is.na(df$umap_2),
]

subtype_order <- c(
  "Normal adjacent tissue",
  "AIS",
  "MIA",
  "IA",
  "COPD",
  "Granulomatosis"
)

subtype_order <- subtype_order[subtype_order %in% unique(df$subtype)]

df$subtype <- factor(df$subtype, levels = subtype_order)

label_df <- df %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

show_subtype_label <- FALSE

p <- ggplot(df, aes(x = umap_1, y = umap_2, colour = subtype)) +
  ggrastr::geom_point_rast(
    size = 0.6,
    stroke = 0,
    shape = 16,
    alpha = 0.9,
    raster.dpi = 300
  )

if (show_subtype_label) {
  p <- p +
    geom_point(
      data = label_df,
      aes(x = umap_1, y = umap_2),
      color = "black",
      fill = "white",
      size = 7,
      alpha = 0.6,
      shape = 21,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_df,
      aes(x = umap_1, y = umap_2, label = subtype),
      colour = "black",
      size = 3.9,
      inherit.aes = FALSE
    )
}

p <- p +
  scale_colour_manual(
    values = subtype_cols,
    breaks = subtype_order,
    name = "Subtype"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 5, alpha = 1)
    )
  ) +
  coord_equal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.key.size = grid::unit(5, "mm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    
    panel.border = element_blank(),
    axis.line = element_blank()
  )

print(p)




