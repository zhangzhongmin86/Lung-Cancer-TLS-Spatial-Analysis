# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)




BiocManager::install("glmGamPoi")



library(future)

plan(multisession, workers = 2)
options(future.globals.maxSize = 8 * 1024^3)
plan(sequential)
####卵巢癌####
BiocManager::install("ensembldb")
BiocManager::install("AnnotationHub")
library(AnnotationHub)
library(BiocFileCache)
library(dbplyr)
library(Seurat)
library(ensembldb)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)


library(Seurat)
library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
library(patchwork)

options(future.globals.maxSize = 2 * 1024^3)

pbmc = Load10X_Spatial(data.dir = "/data/beifen/zhongmin/泛癌空转/卵巢癌/10X/Human Ovarian Cancer, 11 mm Capture Area (FFPE)/原始数据",
                       filename = "CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.h5",
                       assay = "Spatial"
)

plot1 <- VlnPlot(pbmc, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(pbmc, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


pbmc <- SCTransform(pbmc, assay = "Spatial", verbose = FALSE)


SpatialFeaturePlot(pbmc, features = c("CD68", "COL1A2","PECAM1","MMP11","PGF"))

library(ggplot2)
plot <- SpatialFeaturePlot(pbmc, features = c("COL4A1", "COL6A1","LAMA4","ITGA1","ITGB1")) + theme(legend.text = element_text(size = 0),
                                                                                                   legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))

print(plot)


library(Seurat)
library(ggplot2)
library(patchwork)

# 1. 定义你的基因集 (使用你最上方提供的50个基因)
#module_3
gene_list <- c(
  "SPARCL1", "MGP", "CALD1", "IGFBP7", "TIMP3", "A2M", "CAV1", "SPARC", 
  "GNG11", "NFIB", "IGFBP4", "HSPG2", "TM4SF1", "COL4A1", "EPAS1", "LDB2", 
  "COL4A2", "ARHGAP29", "CRIP2", "RAMP2", "CNN3", "CAV2", "EMP2", "EGFL7", 
  "SLC9A3R2", "ENG", "S100A16", "TJP1", "EMP1", "CLEC14A", "NNMT", "TCF4", 
  "PLS3", "CALCRL", "HYAL2", "BCAM", "ESAM", "COX7A1", "JAM2", "FSTL1", 
  "ACVRL1", "LMCD1", "DLC1", "LAMA4", "WWTR1", "APP", "PALMD", "CD34", 
  "PRSS23", "PTPRG"
)

#module_10
gene_list <- c(
  "AIF1", "CD14", "LYZ", "TYROBP", "MS4A6A", "FCER1G", "C5AR1", "MNDA",
  "CLEC7A", "SPI1", "FCGR2A", "CD163", "LST1", "MAFB", "NCF2", "CD68",
  "CSTA", "PILRA", "S100A9", "FPR1", "MS4A7", "CFD", "SLC7A7", "FCGR1A",
  "HCK", "SLC11A1", "CYBB", "CSF1R", "CPVL", "TGFBI", "PLAUR", "CTSS",
  "CST3", "FTL", "HNMT", "S100A8", "TNFSF13B", "PSAP", "MS4A4A", "CSF3R",
  "RAB31", "CD86", "LILRB4", "MPEG1", "CTSB", "TREM1", "IL1B", "FCN1",
  "TLR2", "LRP1"
)

# 2. 检查基因是否在数据中存在 (对应 Python 的检查逻辑)
# 获取当前活跃 Assay (SCT) 中的所有基因名
available_genes <- rownames(pbmc) 
genes_in_data <- intersect(gene_list, available_genes)
genes_not_found <- setdiff(gene_list, available_genes)

if (length(genes_not_found) > 0) {
  cat("以下基因未在数据中找到:\n", paste(genes_not_found, collapse = ", "), "\n\n")
}
cat(sprintf("将使用 %d 个基因进行打分\n", length(genes_in_data)))

# 3. 计算模块打分 (对应 sc.tl.score_genes)
# 注意：features 参数必须是一个 list
pbmc <- AddModuleScore(
  object = pbmc,
  features = list(genes_in_data),
  assay = "SCT",          # 明确使用 SCTransform 后的数据
  name = "Module_Score"   # 设定打分列的名称前缀
)

# ⚠️ 关键提示：Seurat 的 AddModuleScore 会自动在你的 name 后面加上数字序列 (如 1, 2...)
# 因为我们只传入了一个 list，所以结果列名必定是 "Module_Score1"

# 查看打分分布 (对应 print(adata.obs['module10_cell_score'].describe()))
summary(pbmc@meta.data$Module_Score1)

# 4. 可视化：将打分映射到空间位置上 (对应 sc.pl.umap，但是使用空间专属函数)
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1", 
  min.cutoff = 0,
  max.cutoff = 0.7,
  pt.size.factor = 1.6, # 可根据切片点阵的密集程度微调点的大小
  alpha = c(0.1, 1)     # 设定透明度范围，低分透明，高分不透明，有助于看清底层组织结构
) + 
  theme(legend.position = "right") +
  ggtitle("Module_10 Score") +
  scale_fill_viridis_c(option = "C") # 选配：使用类似 Scanpy 常用的高对比度色带 (如 magma)

# 打印图像
print(spatial_score_plot)


library(Seurat)
library(ggplot2)

# 确认模块分数列存在
score_col <- "Module_Score1"
if (!score_col %in% colnames(pbmc@meta.data)) {
  stop(paste0("在 pbmc@meta.data 中没有找到列: ", score_col))
}

# 提取分数
score_vec <- pbmc@meta.data[[score_col]]

# 检查是否全是 NA
if (all(is.na(score_vec))) {
  stop(paste0(score_col, " 全部都是 NA，无法绘图。"))
}

# 计算 1% 和 99% 分位数
q01 <- as.numeric(quantile(score_vec, probs = 0.01, na.rm = TRUE))
q99 <- as.numeric(quantile(score_vec, probs = 0.99, na.rm = TRUE))

cat("Module_Score1 的 1% 分位数:", q01, "\n")
cat("Module_Score1 的 99% 分位数:", q99, "\n")

# 如果分布过窄，避免上下限相同
if (q01 == q99) {
  q01 <- min(score_vec, na.rm = TRUE)
  q99 <- max(score_vec, na.rm = TRUE)
  cat("1% 和 99% 分位数相同，改用最小值和最大值作为上下限。\n")
  cat("新的下限:", q01, "\n")
  cat("新的上限:", q99, "\n")
}

# 将分数截断到 q01 ~ q99，增强可视化对比
pbmc$Module_Score1_clip <- pmax(pmin(score_vec, q99), q01)

# 保持你原来的可视化风格：HE 背景 + SpatialFeaturePlot + viridis 配色
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1_clip",
  min.cutoff = q01,
  max.cutoff = q99,
  pt.size.factor = 1.6,
  alpha = c(0.1, 1)
) +
  theme(legend.position = "right") +
  ggtitle("Module_10 Score") +
  scale_fill_viridis_c(option = "C")

print(spatial_score_plot)


library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

genes_to_plot <- c("COL4A1", "ITGA1")
genes_to_plot <- c("CD68", "LGALS9","TREM2" ,"APOE","TYROBP","SPP1")


# 确保基因存在
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(pbmc)]
print(genes_to_plot)

if (length(genes_to_plot) == 0) {
  stop("这些基因都不在 pbmc 中。")
}

# 每个基因单独作图
plot_list <- lapply(genes_to_plot, function(gene) {
  expr <- FetchData(pbmc, vars = gene)[, 1]
  
  # 去掉 NA
  expr_non_na <- expr[!is.na(expr)]
  
  # 用 1% 和 99% 分位数增强对比
  q_low <- as.numeric(quantile(expr_non_na, probs = 0.01, na.rm = TRUE))
  q_high <- as.numeric(quantile(expr_non_na, probs = 0.99, na.rm = TRUE))
  
  # 如果表达几乎没变化，防止报错
  if (q_low == q_high) {
    q_low <- min(expr_non_na, na.rm = TRUE)
    q_high <- max(expr_non_na, na.rm = TRUE)
  }
  
  p <- SpatialFeaturePlot(
    object = pbmc,
    features = gene,
    min.cutoff = q_low,
    max.cutoff = q_high,
    pt.size.factor = 1.8,
    image.alpha = 0.15,
    crop = TRUE
  ) +
    scale_fill_gradientn(
      colours = c("#313695", "#4575B4", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
      values = rescale(c(
        q_low,
        q_low + (q_high - q_low) * 0.20,
        q_low + (q_high - q_low) * 0.45,
        q_low + (q_high - q_low) * 0.70,
        q_low + (q_high - q_low) * 0.90,
        q_high
      )),
      limits = c(q_low, q_high),
      oob = squish,
      name = gene
    ) +
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(2.2, "cm"),
        barwidth = unit(0.35, "cm"),
        ticks = TRUE
      )
    ) +
    labs(title = NULL) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 7),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )
  
  return(p)
})

# 强制每行 2 个基因，并让子图更紧凑
plot <- wrap_plots(plot_list, ncol = 2, byrow = TRUE)

print(plot)





#####前列腺癌_验证module11#####


library(Seurat)
library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
library(patchwork)




pbmc = Load10X_Spatial(data.dir = "/data/beifen/zhongmin/泛癌空转/前列腺癌/10X/Human Prostate Cancer, Acinar Cell Carcinoma (FFPE)/原始数据",
                       filename = "Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma_filtered_feature_bc_matrix.h5",
                       assay = "Spatial"
)


# pbmc = readRDS("/data/beifen/zhongmin/泛癌空转/前列腺癌/STDS0000047/原始数据/Human-Prostate-Cancer_Adjacent_Normal_Section_10xvisium.rds")



plot1 <- VlnPlot(pbmc, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(pbmc, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


pbmc <- SCTransform(pbmc, assay = "Spatial", verbose = FALSE)


# SpatialFeaturePlot(pbmc, features = c("CD68", "COL1A2","PECAM1","MMP11","PGF"))
# 
library(ggplot2)
plot <- SpatialFeaturePlot(pbmc, features = c("CLEC2C","KLRB1")) + theme(legend.text = element_text(size = 0),
                                                                         legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))

print(plot)


library(Seurat)
library(ggplot2)
library(patchwork)

# 1. 定义你的基因集 (使用你最上方提供的50个基因)
#module_11
gene_list <- c(
  "GZMA", "NKG7", "CD8A", "CTSW", "CD7", "CD2", "CD3E", "GZMM",
  "CST7", "PRF1", "GZMH", "KLRD1", "LCK", "GZMK", "CD3D", "IL32",
  "CD247", "CD96", "KLRB1", "CD8B", "XCL2", "GZMB", "SAMD3", "CD3G",
  "GNLY", "HOPX", "SH2D1A", "PYHIN1", "ZAP70", "IFNG", "IL2RB", "MATK",
  "CD69", "SKAP1", "SH2D2A", "STAT4", "HCST", "XCL1", "SPOCK2", "KLRC1",
  "LINC00861", "LAT", "GPR171", "RUNX3", "SLA2", "ITK", "TIGIT", "BCL11B",
  "SCML4", "ACAP1"
)

# 2. 检查基因是否在数据中存在 (对应 Python 的检查逻辑)
# 获取当前活跃 Assay (SCT) 中的所有基因名
available_genes <- rownames(pbmc) 
genes_in_data <- intersect(gene_list, available_genes)
genes_not_found <- setdiff(gene_list, available_genes)

if (length(genes_not_found) > 0) {
  cat("以下基因未在数据中找到:\n", paste(genes_not_found, collapse = ", "), "\n\n")
}
cat(sprintf("将使用 %d 个基因进行打分\n", length(genes_in_data)))

# 3. 计算模块打分 (对应 sc.tl.score_genes)
# 注意：features 参数必须是一个 list
pbmc <- AddModuleScore(
  object = pbmc,
  features = list(genes_in_data),
  assay = "SCT",          # 明确使用 SCTransform 后的数据
  name = "Module_Score"   # 设定打分列的名称前缀
)

# ⚠️ 关键提示：Seurat 的 AddModuleScore 会自动在你的 name 后面加上数字序列 (如 1, 2...)
# 因为我们只传入了一个 list，所以结果列名必定是 "Module_Score1"

# 查看打分分布 (对应 print(adata.obs['module10_cell_score'].describe()))
summary(pbmc@meta.data$Module_Score1)

# 4. 可视化：将打分映射到空间位置上 (对应 sc.pl.umap，但是使用空间专属函数)
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1", 
  pt.size.factor = 1.6, # 可根据切片点阵的密集程度微调点的大小
  alpha = c(0.1, 1)     # 设定透明度范围，低分透明，高分不透明，有助于看清底层组织结构
) + 
  theme(legend.position = "right") +
  ggtitle("Gene Set Module Score") +
  scale_fill_viridis_c(option = "C") # 选配：使用类似 Scanpy 常用的高对比度色带 (如 magma)

# 打印图像
print(spatial_score_plot)

# 4. 可视化：将打分映射到空间位置上 (对应 sc.pl.umap，但是使用空间专属函数)
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1", 
  min.cutoff = 0,
  max.cutoff = 0.3,
  pt.size.factor = 1.6, # 可根据切片点阵的密集程度微调点的大小
  alpha = c(0.1, 1)     # 设定透明度范围，低分透明，高分不透明，有助于看清底层组织结构
) + 
  theme(legend.position = "right") +
  ggtitle("Module_11 Score") +
  scale_fill_viridis_c(option = "C") # 选配：使用类似 Scanpy 常用的高对比度色带 (如 magma)

# 打印图像
print(spatial_score_plot)



library(Seurat)
library(ggplot2)

# 确认模块分数列存在
score_col <- "Module_Score1"
if (!score_col %in% colnames(pbmc@meta.data)) {
  stop(paste0("在 pbmc@meta.data 中没有找到列: ", score_col))
}

# 提取分数
score_vec <- pbmc@meta.data[[score_col]]

# 检查是否全是 NA
if (all(is.na(score_vec))) {
  stop(paste0(score_col, " 全部都是 NA，无法绘图。"))
}

# 计算 1% 和 99% 分位数
q01 <- as.numeric(quantile(score_vec, probs = 0.01, na.rm = TRUE))
q99 <- as.numeric(quantile(score_vec, probs = 0.99, na.rm = TRUE))

cat("Module_Score1 的 1% 分位数:", q01, "\n")
cat("Module_Score1 的 99% 分位数:", q99, "\n")

# 如果分布过窄，避免上下限相同
if (q01 == q99) {
  q01 <- min(score_vec, na.rm = TRUE)
  q99 <- max(score_vec, na.rm = TRUE)
  cat("1% 和 99% 分位数相同，改用最小值和最大值作为上下限。\n")
  cat("新的下限:", q01, "\n")
  cat("新的上限:", q99, "\n")
}

# 将分数截断到 q01 ~ q99，增强可视化对比
pbmc$Module_Score1_clip <- pmax(pmin(score_vec, q99), q01)

# 保持你原来的可视化风格：HE 背景 + SpatialFeaturePlot + viridis 配色
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1_clip",
  min.cutoff = q01,
  max.cutoff = q99,
  pt.size.factor = 1.6,
  alpha = c(0.1, 1)
) +
  theme(legend.position = "right") +
  ggtitle("Module_9 Score") +
  scale_fill_viridis_c(option = "C")

print(spatial_score_plot)


library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)


genes_to_plot <- c("CLEC2D","CD3D")

# 确保基因存在
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(pbmc)]
print(genes_to_plot)

if (length(genes_to_plot) == 0) {
  stop("这些基因都不在 pbmc 中。")
}

# 每个基因单独作图
plot_list <- lapply(genes_to_plot, function(gene) {
  expr <- FetchData(pbmc, vars = gene)[, 1]
  
  # 去掉 NA
  expr_non_na <- expr[!is.na(expr)]
  
  # 用 1% 和 99% 分位数增强对比
  q_low <- as.numeric(quantile(expr_non_na, probs = 0.01, na.rm = TRUE))
  q_high <- as.numeric(quantile(expr_non_na, probs = 0.99, na.rm = TRUE))
  
  # 如果表达几乎没变化，防止报错
  if (q_low == q_high) {
    q_low <- min(expr_non_na, na.rm = TRUE)
    q_high <- max(expr_non_na, na.rm = TRUE)
  }
  
  p <- SpatialFeaturePlot(
    object = pbmc,
    features = gene,
    min.cutoff = q_low,
    max.cutoff = q_high,
    pt.size.factor = 1.8,
    image.alpha = 0.15,
    crop = TRUE
  ) +
    scale_fill_gradientn(
      colours = c("#313695", "#4575B4", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
      values = rescale(c(
        q_low,
        q_low + (q_high - q_low) * 0.20,
        q_low + (q_high - q_low) * 0.45,
        q_low + (q_high - q_low) * 0.70,
        q_low + (q_high - q_low) * 0.90,
        q_high
      )),
      limits = c(q_low, q_high),
      oob = squish,
      name = gene
    ) +
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(2.2, "cm"),
        barwidth = unit(0.35, "cm"),
        ticks = TRUE
      )
    ) +
    labs(title = NULL) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 7),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )
  
  return(p)
})

# 强制每行 2 个基因，并让子图更紧凑
plot <- wrap_plots(plot_list, ncol = 2, byrow = TRUE)

print(plot)



#####大肠癌_验证module6#####


library(Seurat)
library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
library(patchwork)




pbmc = Load10X_Spatial(data.dir = "/data/beifen/zhongmin/泛癌空转/大肠癌/10X/Human Intestine Cancer (FPPE)/原始数据",
                       filename = "Visium_FFPE_Human_Intestinal_Cancer_filtered_feature_bc_matrix.h5",
                       assay = "Spatial"
)


# pbmc = readRDS("/data/beifen/zhongmin/泛癌空转/前列腺癌/STDS0000047/原始数据/Human-Prostate-Cancer_Adjacent_Normal_Section_10xvisium.rds")



plot1 <- VlnPlot(pbmc, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(pbmc, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


pbmc <- SCTransform(pbmc, assay = "Spatial", verbose = FALSE)


# SpatialFeaturePlot(pbmc, features = c("CD68", "COL1A2","PECAM1","MMP11","PGF"))
# 
# library(ggplot2)
# plot <- SpatialFeaturePlot(pbmc, features = c("CD68", "COL1A2","PECAM1","MMP11","PGF")) + theme(legend.text = element_text(size = 0),
#                                                                                                 legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
# 
# print(plot)


library(Seurat)
library(ggplot2)
library(patchwork)

# 1. 定义你的基因集 (使用你最上方提供的50个基因)
#module_6
gene_list <- c(
  "CAV1", "SPARCL1", "SPARC", "TCF4", "RAMP2", "CLIC4", "CAV2", "MGP",
  "CNN3", "CALD1", "COL4A1", "A2M", "TIMP3", "HSPG2", "TM4SF1", "CLEC14A",
  "COL4A2", "GAS6", "LAMC1", "NFIB", "S100A13", "SLC9A3R2", "IGFBP7", "IGFBP4",
  "GNG11", "NNMT", "LAMA5", "JAM2", "PDLIM1", "SERPINH1", "ID1", "LDB2",
  "PRDX4", "LMCD1", "ITGA6", "KCNN3", "TJP1", "EMP2", "ADIRF", "LAMA4",
  "EGFL7", "CALCRL", "CTTN", "NES", "HYAL2", "RHOJ", "RAB13", "ARHGAP29",
  "ID3", "COL15A1"
)

# 2. 检查基因是否在数据中存在 (对应 Python 的检查逻辑)
# 获取当前活跃 Assay (SCT) 中的所有基因名
available_genes <- rownames(pbmc) 
genes_in_data <- intersect(gene_list, available_genes)
genes_not_found <- setdiff(gene_list, available_genes)

if (length(genes_not_found) > 0) {
  cat("以下基因未在数据中找到:\n", paste(genes_not_found, collapse = ", "), "\n\n")
}
cat(sprintf("将使用 %d 个基因进行打分\n", length(genes_in_data)))

# 3. 计算模块打分 (对应 sc.tl.score_genes)
# 注意：features 参数必须是一个 list
pbmc <- AddModuleScore(
  object = pbmc,
  features = list(genes_in_data),
  assay = "SCT",          # 明确使用 SCTransform 后的数据
  name = "Module_Score"   # 设定打分列的名称前缀
)

# ⚠️ 关键提示：Seurat 的 AddModuleScore 会自动在你的 name 后面加上数字序列 (如 1, 2...)
# 因为我们只传入了一个 list，所以结果列名必定是 "Module_Score1"

# 查看打分分布 (对应 print(adata.obs['module10_cell_score'].describe()))
summary(pbmc@meta.data$Module_Score1)

# 4. 可视化：将打分映射到空间位置上 (对应 sc.pl.umap，但是使用空间专属函数)
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1", 
  pt.size.factor = 1.6, # 可根据切片点阵的密集程度微调点的大小
  alpha = c(0.1, 1)     # 设定透明度范围，低分透明，高分不透明，有助于看清底层组织结构
) + 
  theme(legend.position = "right") +
  ggtitle("Gene Set Module Score") +
  scale_fill_viridis_c(option = "C") # 选配：使用类似 Scanpy 常用的高对比度色带 (如 magma)

# 打印图像
print(spatial_score_plot)


# 4. 可视化：将打分映射到空间位置上
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1", 
  min.cutoff = 0,
  max.cutoff = 0.8,
  pt.size.factor = 1.6, # 可根据切片点阵的密集程度微调点的大小
  alpha = c(0.1, 1)     # 设定透明度范围，低分透明，高分不透明，有助于看清底层组织结构
) + 
  theme(legend.position = "right") +
  ggtitle("Module_6 Score") +
  scale_fill_viridis_c(option = "C") # 选配：使用类似 Scanpy 常用的高对比度色带 (如 magma)

# 打印图像
print(spatial_score_plot)


library(Seurat)
library(ggplot2)

# 确认模块分数列存在
score_col <- "Module_Score1"
if (!score_col %in% colnames(pbmc@meta.data)) {
  stop(paste0("在 pbmc@meta.data 中没有找到列: ", score_col))
}

# 提取分数
score_vec <- pbmc@meta.data[[score_col]]

# 检查是否全是 NA
if (all(is.na(score_vec))) {
  stop(paste0(score_col, " 全部都是 NA，无法绘图。"))
}

# 计算 1% 和 99% 分位数
q01 <- as.numeric(quantile(score_vec, probs = 0.01, na.rm = TRUE))
q99 <- as.numeric(quantile(score_vec, probs = 0.99, na.rm = TRUE))

cat("Module_Score1 的 1% 分位数:", q01, "\n")
cat("Module_Score1 的 99% 分位数:", q99, "\n")

# 如果分布过窄，避免上下限相同
if (q01 == q99) {
  q01 <- min(score_vec, na.rm = TRUE)
  q99 <- max(score_vec, na.rm = TRUE)
  cat("1% 和 99% 分位数相同，改用最小值和最大值作为上下限。\n")
  cat("新的下限:", q01, "\n")
  cat("新的上限:", q99, "\n")
}

# 将分数截断到 q01 ~ q99，增强可视化对比
pbmc$Module_Score1_clip <- pmax(pmin(score_vec, q99), q01)

# 保持你原来的可视化风格：HE 背景 + SpatialFeaturePlot + viridis 配色
spatial_score_plot <- SpatialFeaturePlot(
  object = pbmc,
  features = "Module_Score1_clip",
  min.cutoff = q01,
  max.cutoff = q99,
  pt.size.factor = 1.6,
  alpha = c(0.1, 1)
) +
  theme(legend.position = "right") +
  ggtitle("Module_10 Score") +
  scale_fill_viridis_c(option = "C")

print(spatial_score_plot)


library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)


genes_to_plot <- c( "COL1A2","CD74","PECAM1","MZB1")



# 确保基因存在
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(pbmc)]
print(genes_to_plot)

if (length(genes_to_plot) == 0) {
  stop("这些基因都不在 pbmc 中。")
}

# 每个基因单独作图
plot_list <- lapply(genes_to_plot, function(gene) {
  expr <- FetchData(pbmc, vars = gene)[, 1]
  
  # 去掉 NA
  expr_non_na <- expr[!is.na(expr)]
  
  # 用 1% 和 99% 分位数增强对比
  q_low <- as.numeric(quantile(expr_non_na, probs = 0.01, na.rm = TRUE))
  q_high <- as.numeric(quantile(expr_non_na, probs = 0.99, na.rm = TRUE))
  
  # 如果表达几乎没变化，防止报错
  if (q_low == q_high) {
    q_low <- min(expr_non_na, na.rm = TRUE)
    q_high <- max(expr_non_na, na.rm = TRUE)
  }
  
  p <- SpatialFeaturePlot(
    object = pbmc,
    features = gene,
    min.cutoff = q_low,
    max.cutoff = q_high,
    pt.size.factor = 1.8,
    image.alpha = 0.15,
    crop = TRUE
  ) +
    scale_fill_gradientn(
      colours = c("#313695", "#4575B4", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
      values = rescale(c(
        q_low,
        q_low + (q_high - q_low) * 0.20,
        q_low + (q_high - q_low) * 0.45,
        q_low + (q_high - q_low) * 0.70,
        q_low + (q_high - q_low) * 0.90,
        q_high
      )),
      limits = c(q_low, q_high),
      oob = squish,
      name = gene
    ) +
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(2.2, "cm"),
        barwidth = unit(0.35, "cm"),
        ticks = TRUE
      )
    ) +
    labs(title = NULL) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 7),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )
  
  return(p)
})

# 强制每行 2 个基因，并让子图更紧凑
plot <- wrap_plots(plot_list, ncol = 2, byrow = TRUE)

print(plot)


