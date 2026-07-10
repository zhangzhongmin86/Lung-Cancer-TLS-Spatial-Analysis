# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)




####NMF非负矩阵分解######
library(registry)
library(registry)
library(cluster)

library(glue)
library(dplyr)
library(NMF)
library(reshape2)
library(ComplexHeatmap)


obs_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "数据", "cellchat_communication_matrix.csv")#group列有NA值
info <- read.csv(obs_file, row.names = 1, check.names = TRUE)
table(info$tissue_organ)
head(info[1:5,1:5])

rows_to_keep <- !grepl("HLA-A|HLA-B|HLA-C|HLA-D|HLA-E|HLA-F", rownames(info))
info <- info[rows_to_keep, ]

sample.info <- as.data.frame(readr::read_csv(file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "cellchat_按患者分.csv"), show_col_types = FALSE))
head(sample.info[1:5,1:5])

cat("sample.info中有", length(unique(sample.info$sample)), "个唯一sample\n")
table(sample.info$tissue_organ)


ratio = as.data.frame(t(info))
ratio$sample = row.names(ratio)

# 假设 ratio$sample 已经赋值了 row.names(ratio)
ratio$sample <- as.character(ratio$sample)

# 去掉 sample 列中，第一位是 X 且第二位是数字(0-9) 的前缀 X
ratio$sample <- sub("^X([0-9])", "\\1", ratio$sample)
ratio$sample <- gsub("\\.", "-", ratio$sample)
# 1. 将特定样本名中的连字符(-)改为点号(.)
ratio <- ratio %>%
  mutate(sample = case_when(
    sample == "NP44-NB_v1-1" ~ "NP44-NB_v1.1",
    sample == "NP44-NB_v1-0" ~ "NP44-NB_v1.0",
    sample == "HHD-1" ~ "HHD.1",
    TRUE ~ sample  # 其他样本保持不变
  ))

# 2. 将以"SSC-"开头的样本名中的连字符(-)改为空格
ratio$sample <- ifelse(startsWith(ratio$sample, "SSC-"),
                       sub("^SSC-", "SSC ", ratio$sample),
                       ratio$sample)

# 3. 将包含"Normal-adjacent-tissue"的样本名中的连字符(-)改为空格
ratio$sample <- ifelse(grepl("Normal-adjacent-tissue", ratio$sample, fixed = TRUE),
                       gsub("-", " ", ratio$sample, fixed = TRUE),
                       ratio$sample)
# 行名也做同样处理（如果行名用作ID很重要的话！）
rownames(ratio) <- ratio$sample


sample.info <- sample.info[sample.info$sample %in% ratio$sample,]


response.meta <- sample.info
head(response.meta[1:5,1:5])

merge.version <- merge(ratio, response.meta, by = "sample", all.x = TRUE)

cat("ratio中有", length(unique(ratio$sample)), "个唯一sample\n")
cat("成功对应上的sample数量为", length(intersect(ratio$sample, sample.info$sample)), "\n")


not_in_sampleinfo <- setdiff(ratio$sample, sample.info$sample)
cat("ratio中有", length(not_in_sampleinfo), "个sample在sample.info中没有找到:\n")
print(not_in_sampleinfo[1:20])  # 只显示前10个，防止太长


ratio <- ratio[merge.version$sample,]
ratio$sample <- NULL
ratio[is.na(ratio)] <- 0

# 3. 清理全0行与全0列
info_mat <- as.matrix(ratio)
row_zeros <- rowSums(info_mat) == 0
col_zeros <- colSums(info_mat) == 0
info_mat2 <- info_mat[!row_zeros, !col_zeros]
info_clean <- as.data.frame(info_mat2)
cat("过滤全零行列后，维度为：", paste(dim(info_clean), collapse = " x "), "\n")

# 4. 保留“在3个上样本有信号”的通讯对
nonzero_count_per_col <- colSums(info_clean != 0)
cols_10plus_mat <- info_clean[, nonzero_count_per_col >= 3]
cat("筛选后列数（3行及以上有信号）：", ncol(cols_10plus_mat), "\n")


# 5. 可选：进一步feature reduction，保留方差最高前10000个通讯对
var_col <- apply(cols_10plus_mat, 2, var)
topN <- 10000
select_col_idx <- order(var_col, decreasing = TRUE)[1:min(topN, ncol(cols_10plus_mat))]
cols_var_top <- cols_10plus_mat[, select_col_idx]
cat("最终用于NMF的疾病数量（top方差）：", ncol(cols_var_top), "\n")

# output_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "preprocessed_communication_matrix_cols_var_top_10000.csv")
# write.csv(cols_var_top, file = output_file, row.names = TRUE)

preprocessed_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "preprocessed_communication_matrix_cols_var_top_10000.csv")
cols_var_top <- read.csv(preprocessed_file, row.names = 1, check.names = FALSE)

rows_var_top_log <- apply(cols_var_top, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
rows_var_top_log <- as.data.frame(rows_var_top_log)
rows_var_top_log <- t(rows_var_top_log)
head(rows_var_top_log[1:5,1:5])

# output_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "preprocessed_communication_matrix_10000.csv")
# write.csv(rows_var_top_log, file = output_file, row.names = TRUE)

preprocessed_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "preprocessed_communication_matrix_10000.csv")
# 读取时注意保持行名
rows_var_top_log <- read.csv(preprocessed_file, row.names = 1, check.names = FALSE)
# 转换为矩阵格式（NMF需要矩阵输入）
rows_var_top_log <- as.matrix(rows_var_top_log)
head(rows_var_top_log[1:5,1:5])
dim(rows_var_top_log)
###

set.seed(100)
ranks <- 2:20
estim.coad <- nmf(rows_var_top_log, ranks, nrun = 2, method = "lee")
plot(estim.coad)
saveRDS(estim.coad, file = file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat", "结果", "nmf_result_cellchat.rds"))

estim.coad = readRDS(file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat", "结果", "nmf_result_cellchat.rds"))
plot(estim.coad)

# 8. 选定一个rank，做最终分解
best_rank <- 14 # 可根据上面结果调整
final.nmf <- nmf(rows_var_top_log, best_rank, nrun = 100, method = "lee")
saveRDS(final.nmf, file = file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "nmf_result_cellchat_best14_2_1000.rds"))

###绘图####


estim.coad <- readRDS(file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "结果", "去除HLA", "nmf_result_cellchat_去除HLA.rds"))
plot(estim.coad)
dim(estim.coad)



nmf.rank5 = readRDS(file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "nmf_result_cellchat_best14_1000.rds"))

plot(nmf.rank5)
dim(nmf.rank5)

index <- extractFeatures(nmf.rank5, "max")  # change the order of the index

new.index <- list()
new.index[[1]] <- index[[1]]
new.index[[2]] <- index[[2]]
new.index[[3]] <- index[[3]]
new.index[[4]] <- index[[4]]
new.index[[5]] <- index[[5]]
new.index[[6]] <- index[[6]]
new.index[[7]] <- index[[7]]
new.index[[8]] <- index[[8]]
new.index[[9]] <- index[[9]]
new.index[[10]] <- index[[10]]
new.index[[11]] <- index[[11]]
new.index[[12]] <- index[[12]]
new.index[[13]] <- index[[13]]
new.index[[14]] <- index[[14]]
sapply(index, length)


sig.order <- unlist(new.index)

NMF.Exp.rank5 <- rows_var_top_log[sig.order,]
# NMF.Exp.rank5 <- na.omit(NMF.Exp.rank5)
NMF.Exp.rank5 <- na.omit(NMF.Exp.rank5)

group <- predict(nmf.rank5)  # adjust the position of the module

new.group <- c()
for (each in group) {
  if (each %in% c("1")) { new.group <- c(new.group, "1") }
  if (each %in% c("2")) { new.group <- c(new.group, "2") }
  if (each %in% c("3")) { new.group <- c(new.group, "3") }
  if (each %in% c("4")) { new.group <- c(new.group, "4") }
  if (each %in% c("5")) { new.group <- c(new.group, "5") }
  if (each %in% c("6")) { new.group <- c(new.group, "6") }
  if (each %in% c("7")) { new.group <- c(new.group, "7") }
  if (each %in% c("8")) { new.group <- c(new.group, "8") }
  if (each %in% c("9")) { new.group <- c(new.group, "9") }
  if (each %in% c("10")) { new.group <- c(new.group, "10") }
  if (each %in% c("11")) { new.group <- c(new.group, "11") }
  if (each %in% c("12")) { new.group <- c(new.group, "12") }
  if (each %in% c("13")) { new.group <- c(new.group, "13") }
  if (each %in% c("14")) { new.group <- c(new.group, "14") }
}

new.group <- factor(new.group, levels = c("1", "2", "3", "4","5","6","7","8","9","10","11","12","13","14"))#,"14"
new.group <- factor(new.group, levels = c("1", "2", "3", "4","5","6","7","8","9","10","11","12","13"))#,"14"

z_ratio <-scale(cols_var_top)/4
z_ratio <-as.data.frame(z_ratio)
z_ratio <-t(z_ratio)
plot_matrix <- z_ratio[sig.order,]
plot_matrix <- na.omit(plot_matrix)


plot_matrix <- z_ratio[sig.order, , drop = FALSE]

# 选项1：把 NA 置 0（或置为行均值/中位数）
plot_matrix[is.na(plot_matrix)] <- 0

info.matrix <-as.data.frame(t(NMF.Exp.rank5))
info.matrix$sample<-rownames(info.matrix)
info.matrix$group<-new.group
info.matrix <-merge(info.matrix, response.meta, by ="sample", all.x =TRUE)

# output_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "nmf结果14.csv")
# write.csv(info.matrix, file = output_file, row.names = TRUE)

dim(info.matrix)

gene.group <- c()
for (each in rownames(NMF.Exp.rank5)) {
  if (each %in% rownames(rows_var_top_log)[new.index[[1]]]) {
    gene.group <- c(gene.group, "module1")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[2]]]) {
    gene.group <- c(gene.group, "module2")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[3]]]) {
    gene.group <- c(gene.group, "module3")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[4]]]) {
    gene.group <- c(gene.group, "module4")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[5]]]) {
    gene.group <- c(gene.group, "module5")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[6]]]) {
    gene.group <- c(gene.group, "module6")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[7]]]) {
    gene.group <- c(gene.group, "module7")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[8]]]) {
    gene.group <- c(gene.group, "module8")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[9]]]) {
    gene.group <- c(gene.group, "module9")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[10]]]) {
    gene.group <- c(gene.group, "module10")
  }  else if (each %in% rownames(rows_var_top_log)[new.index[[11]]]) {
    gene.group <- c(gene.group, "module11")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[12]]]) {
    gene.group <- c(gene.group, "module12")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[13]]]) {
    gene.group <- c(gene.group, "module13")
  } else if (each %in% rownames(rows_var_top_log)[new.index[[14]]]) {
    gene.group <- c(gene.group, "module14")
  } 
}


table(gene.group)


table(info.matrix$group.x, useNA = "ifany")
table(info.matrix$group.y, useNA = "ifany")

ha <- HeatmapAnnotation(
  group = factor(info.matrix$group.y, levels = c("Autoimmune diseases", "Inflammation", "Normal","Normal adjacent tissue","Precancerous condition","Tumor","Blood_Tumor","Tumor_metastasis")),
  tissue_organ = factor(info.matrix$tissue_organ, levels = c("adnexa of uterus", "Adrenal", "Bladder", "Blood", "Blood_vessel", "Bone", "Bone_marrow", "Brain", "breast", "Esophagus", 
                                                             "Intestine", "Kidney", "Liver", "Lung", "Lymph_node", "Nose", "Oral_cavity", "Other", "Ovary", "prostate", 
                                                             "Skin", "spleen", "Stomach", "Uterus"
  )),
  cancer = factor(info.matrix$cancer, levels = c("ALC", "atherosclerosis", "atrophic gastritis", "B-ALL", "BE", "BLCA", "BRCA", "CD", "CESC", "CF", "CLL", "COPD", "COVID-19", 
                                                 "CR", "CRC", "cSCC", "CTCL", "ESCC", "FL", "Gastritis", "GBM", "GIM", "HCC", "HLE", "HNSCC", "ICC", "ILD", 
                                                 "KIRC", "KO", "LAM", "LBTC", "LCC", "LCNEC", "light chain disease", "LUAD", "LUAD+LCNEC", "LUSC", "MM", "Normal", "NPC", 
                                                 "NSCLC", "OS", "OSCC", "OV", "PBC", "PC", "PDAC", "PF", "PNA", "PRAD", "PSC", "RA", "sAH", "SCLC", "SS", "STAD", "UC", "UCEC"
  )),
  
  
  module = factor(info.matrix$group.x, levels = c("1", "2", "3", "4","5","6","7","8","9","10","11","12","13","14")),
  col = list(
    group = c("Autoimmune diseases" = "#E5D2DD", "Inflammation" = "#53A85F", "Normal" = "#F1BB72", 
              "Normal adjacent tissue" = "#F3B1A0", "Precancerous condition" = "#D6E7A3", "Tumor" = "#57C3F3", 
              "Blood_Tumor" = "#476D87", "Tumor_metastasis" = "#E95C59"),
    tissue_organ = c("adnexa of uterus" = "#E59CC4", "Adrenal" = "#AB3282", "Bladder" = "#23452F", "Blood" = "#BD956A", 
                     "Blood_vessel" = "#8C549C", "Bone" = "#585658", "Bone_marrow" = "#9FA3A8", "Brain" = "#E0D4CA", 
                     "breast" = "#5F3D69", "Esophagus" = "#C5DEBA", "Intestine" = "#58A4C3", "Kidney" = "#E4C755", 
                     "Liver" = "#F7F398", "Lung" = "#AA9A59", "Lymph_node" = "#E63863", "Nose" = "#E39A35", 
                     "Oral_cavity" = "#C1E6F3", "Other" = "#6778AE", "Ovary" = "#91D0BE", "prostate" = "#B53E2B", 
                     "Skin" = "#712820", "spleen" = "#DCC1DD", "Stomach" = "#CCE0F5", "Uterus" = "#CCC9E6"),
    cancer = c("ALC" = "#E8F3F1", "atherosclerosis" = "#F2F1EB", "atrophic gastritis" = "#EEE7DA", "B-ALL" = "#AFC8AD", 
               "BE" = "#88AB8E", "BLCA" = "#E97777", "BRCA" = "#FF9F9F", "CD" = "#FCDDB0", "CESC" = "#E5D2DD", 
               "CF" = "#53A85F", "CLL" = "#F1BB72", "COPD" = "#F3B1A0", "COVID-19" = "#D6E7A3", "CR" = "#57C3F3", 
               "CRC" = "#476D87", "cSCC" = "#E95C59", "CTCL" = "#E59CC4", "ESCC" = "#AB3282", "FL" = "#23452F", 
               "Gastritis" = "#BD956A", "GBM" = "#8C549C", "GIM" = "#585658", "HCC" = "#9FA3A8", "HLE" = "#E0D4CA", 
               "HNSCC" = "#5F3D69", "ICC" = "#C5DEBA", "ILD" = "#58A4C3", "KIRC" = "#E4C755", "KO" = "#F7F398", 
               "LAM" = "#AA9A59", "LBTC" = "#E63863", "LCC" = "#E39A35", "LCNEC" = "#C1E6F3", "light chain disease" = "#6778AE", 
               "LUAD" = "#91D0BE", "LUAD+LCNEC" = "#B53E2B", "LUSC" = "#712820", "MM" = "#DCC1DD", "Normal" = "#CCE0F5", 
               "NPC" = "#CCC9E6", "NSCLC" = "#625D9E", "OS" = "#68A180", "OSCC" = "#3A6963", "OV" = "#968175", 
               "PBC" = "#E8F3F1", "PC" = "#F2F1EB", "PDAC" = "#EEE7DA", "PF" = "#AFC8AD", "PNA" = "#88AB8E", "PRAD" = "#E97777", 
               "PSC" = "#FF9F9F", "RA" = "#FCDDB0", "sAH" = "#E5D2DD", "SCLC" = "#53A85F", "SS" = "#F1BB72", "STAD" = "#F3B1A0", 
               "UC" = "#D6E7A3", "UCEC" = "#57C3F3"),
    
    module = c("1" = "#E64B35B2", "2" = "#4DBBD5B2", "3" = "#00A087B2", "4" = "#3C5488B2", "5" = "#E5D2DD", "6" ="#53A85F","7" ="#F1BB72","8" = "#57C3F3", 
               "9" = "#E8F3F1", "10" =  "#9FA3A8","11" = "#D6E7A3","12"= "#E95C59","13"= "#AB3282","14"= "#23452F"  )#, "5" = "#F39B7FB2"
  ),
  simple_anno_size = unit(0.5, "cm")
)

library(circlize)

# 创建右侧颜色分组
gene_group_colors <- c(
  "module1" = "#E37D6F", 
  "module2" = "#81C8D9", 
  "module3" = "#129982", 
  "module4" = "#7584A7", 
  "module5" = "#F1B4A3"
)

# 将分组内的每个基因按照上述分组赋值
text_colors <- gene_group_colors[gene.group]

# 添加行注释（右侧），使文本（基因名）被彩色背景框包围
ra <- rowAnnotation(
  gene = anno_text(
    rownames(plot_matrix), 
    gp = gpar(
      fontsize = 10, 
      col = "black", 
      fill = text_colors, 
      box_lwd = 1, 
      box_col = "black"
    ), 
    show_name = FALSE
  ),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.5, "cm")
)

col_fun = colorRamp2(
  breaks = c(0, 0.2, 0.8), # 自定义数据分段点
  # colors = c("#4575b4",  "#fee090", "#d73027")
  colors = c("#4575b4",  '#D6E7A3', '#E95C59')
)


hp <- Heatmap(
  plot_matrix,
  col = col_fun,
  name = "ratio",
  top_annotation = ha,
  # right_annotation = ra,
  row_split = gene.group,
  column_split = new.group,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = order(factor(
    info.matrix$group.y,
    levels = c("Autoimmune diseases","Inflammation","Normal","Normal adjacent tissue",
               "Precancerous condition","Tumor","Tumor_blood","Tumor_metastasis")
  )),
  row_names_gp = grid::gpar(fontsize = 10),
  show_column_names = FALSE,
  show_row_names = FALSE,
  use_raster = FALSE,
  rect_gp = grid::gpar(col = NA, lwd = 0),
  border = NA,
  
  ## 关键：让 module 标签横着放
  row_title_rot = 0,                 # ← 横向
  row_title_side = "left",           # ← 放在左边（可选，默认就是 left）
  row_title_gp = grid::gpar(fontsize = 10)  # ← 大小/样式（可选）
)
draw(hp)

# 使用 ggplot2 将热图保存为文件
library(ggplot2)
heatmap_grob <- grid::grid.grab()
heatmap_ggplot <- cowplot::ggdraw() + cowplot::draw_grob(heatmap_grob)

ggsave(heatmap_ggplot, filename = file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "cell_heatmap14_2025_11_03.png"), width = 15, height = 15, dpi = 600, device = 'png', bg = '#FFFFFF')
ggsave(heatmap_ggplot, filename = file.path(PROJECT_ROOT, "泛癌S", "NMF", "cellchat_sample", "重新提取cellchat结果", "cell_heatmap14_2025_11_03.pdf"), width = 15, height = 15, dpi = 50, device = 'pdf', bg = '#FFFFFF')





####个module种的分组占比######


library(dplyr)
library(ggplot2)

table(info.matrix$group.x, useNA = "ifany")
info.matrix <- info.matrix[!is.na(info.matrix$group.x), ]
filtered <- info.matrix

proportion_data <- filtered %>%
  dplyr::group_by(group.x, group.y) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(group.x) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup()

# 关键修改：反转 group.x 顺序
# 这样在 coord_flip() 后，1 会在最上面，13 会在最下面
proportion_data$group.x <- factor(
  proportion_data$group.x,
  levels = rev(sort(unique(proportion_data$group.x)))
)

# group.y 颜色
group_y_colors <- c(
  "Autoimmune diseases"    = "#F1BB72",
  "Blood_Tumor"            = "#E95C59",
  "Inflammation"           = "#D6E7A3",
  "Normal"                 = "#E5D2DD",
  "Normal adjacent tissue" = "#57C3F3",
  "Precancerous condition" = "#E59CC4",
  "Tumor"                  = "#53A85F",
  "Tumor_metastasis"       = "#F3B1A0",
  "Tumor_blood"            = "#476D87"
)

# 图例和堆叠顺序
group_y_order <- c(
  "Tumor_metastasis",
  "Tumor",
  "Precancerous condition",
  "Normal adjacent tissue",
  "Normal",
  "Inflammation",
  "Blood_Tumor",
  "Autoimmune diseases"
)

proportion_data$group.y <- factor(proportion_data$group.y, levels = group_y_order)

# 标签：仅显示 >1.5% 的占比
proportion_data$label <- ifelse(
  proportion_data$percentage > 5,
  sprintf("%.1f%%", proportion_data$percentage),
  ""
)

p <- ggplot(proportion_data, aes(x = group.x, y = percentage, fill = group.y)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  scale_fill_manual(values = group_y_colors, drop = FALSE) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  coord_flip() +
  labs(
    title = "Percentage of  Categories within Each Module",
    x = "Modules",
    y = "Percentage (%)",
    fill = "Group Y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)





####热图 + 层次聚类——疾病——改R/O比值#####

library(dplyr)
library(pheatmap)

names(info.matrix) <- make.names(names(info.matrix), unique = TRUE)

required_cols <- c("group.x", "cancer", "group.y")
if (!all(required_cols %in% colnames(info.matrix))) {
  stop("info.matrix 中缺少必须的列：group.x / cancer / group.y")
}

filtered2 <- info.matrix %>%
  mutate(
    group.x = as.character(group.x),
    cancer = as.character(cancer),
    group.y = as.character(group.y)
  ) %>%
  filter(!is.na(group.x), !is.na(cancer), !is.na(group.y)) %>%
  filter(group.x != "14") %>%
  filter(group.y != "Normal adjacent tissue") %>%
  filter(cancer != "")

cancer_counts <- table(filtered2$cancer)
small_cancers <- names(cancer_counts)[cancer_counts < 3]

if (length(small_cancers) > 0) {
  message("将被过滤的疾病: ", paste(small_cancers, collapse = ", "))
  filtered2 <- filtered2[!filtered2$cancer %in% small_cancers, , drop = FALSE]
}

obs_mat <- table(filtered2$cancer, filtered2$group.x)
obs_mat <- as.matrix(obs_mat)

colnames(obs_mat) <- paste0("Module_", colnames(obs_mat))

desired_modules <- paste0("Module_", 1:13)
keep_modules <- intersect(desired_modules, colnames(obs_mat))

if (length(keep_modules) == 0) {
  stop("过滤后没有 Module_1 ~ Module_13 的列，请检查 group.x")
}

obs_mat <- obs_mat[, keep_modules, drop = FALSE]
obs_mat <- obs_mat[rowSums(obs_mat) > 0, , drop = FALSE]
obs_mat <- obs_mat[, colSums(obs_mat) > 0, drop = FALSE]

exp_mat <- outer(rowSums(obs_mat), colSums(obs_mat)) / sum(obs_mat)

log2_oe <- log2((obs_mat + 0.5) / (exp_mat + 0.5))

cap_value <- quantile(abs(log2_oe), 0.98, na.rm = TRUE)
cap_value <- max(cap_value, 1)
cap_value <- min(cap_value, 3)

log2_oe_cap <- pmax(pmin(log2_oe, cap_value), -cap_value)

dominant_idx <- max.col(log2_oe, ties.method = "first")
dominant_module <- colnames(log2_oe)[dominant_idx]
dominant_num <- as.integer(sub("Module_", "", dominant_module))

row_names_all <- rownames(log2_oe)

row_order_list <- lapply(sort(unique(dominant_num)), function(m) {
  rn <- row_names_all[dominant_num == m]
  
  if (length(rn) == 1) {
    return(rn)
  }
  
  submat <- log2_oe_cap[rn, , drop = FALSE]
  hc <- hclust(dist(submat), method = "ward.D2")
  rn <- rn[hc$order]
  
  rn <- rn[order(log2_oe[rn, paste0("Module_", m)], decreasing = TRUE)]
  
  rn
})

row_order <- unlist(row_order_list, use.names = FALSE)
plot_mat <- log2_oe_cap[row_order, , drop = FALSE]

group_sizes <- sapply(row_order_list, length)
gaps_row <- cumsum(group_sizes)
gaps_row <- gaps_row[gaps_row < nrow(plot_mat)]

my_palette <- colorRampPalette(c(
  "#476D87",
  "#AFC0CF",
  "#F7F7F7",
  "#F3B8B1",
  "#E95C59"
))(100)

my_breaks <- seq(-cap_value, cap_value, length.out = 101)

legend_ticks <- c(-cap_value, -cap_value / 2, 0, cap_value / 2, cap_value)
legend_labels <- sprintf("%.2f", legend_ticks)

pheatmap(
  plot_mat,
  color = my_palette,
  breaks = my_breaks,
  main = "Cancer-module log2(O/E) heatmap",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_row,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 11,
  angle_col = 90,
  cellwidth = 16,
  cellheight = 12,
  display_numbers = FALSE,
  legend_breaks = legend_ticks,
  legend_labels = legend_labels,
  na_col = "white"
)






