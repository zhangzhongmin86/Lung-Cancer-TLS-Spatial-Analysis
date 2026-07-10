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


obs_file <- file.path(PROJECT_ROOT, "泛癌S", "NMF", "基质免疫样本2", "samples_with_Immune cells_其他另外一种.csv")#group列有NA值
info <- read.csv(obs_file, row.names = 1, check.names = TRUE)
table(info$tissue_organ)

sample.info = info

df <- table(info$sample, info$celltype_3)

ratio <- as.data.frame(df / rowSums(df))
colnames(ratio) <- c("sample", "celltype_3", "Freq")

table(sample.info$group)


pathological_response_level <- c()

library(dplyr)

sample.info <- sample.info %>%
  mutate(pathological_response_level = case_when(
    group == "Autoimmune diseases" ~ "Autoimmune diseases",
    group == "Blood_Tumor" ~ "Blood_Tumor",
    group == "Inflammation" ~ "Inflammation",
    group == "Normal" ~ "Normal",
    group == "Normal adjacent tissue" ~ "Normal adjacent tissue",
    group == "Precancerous condition" ~ "Precancerous condition",
    group == "Tumor" ~ "Tumor",
    group == "Tumor_blood" ~ "Tumor_blood",
    group == "Tumor_metastasis" ~ "Tumor_metastasis"
  ))

# 查看结果
table(sample.info$pathological_response_level)

response.meta = sample.info

response.meta <- response.meta %>% distinct(sample, .keep_all = TRUE)



ratio <- reshape2::dcast(ratio, sample ~ celltype_3, value.var = "Freq")

merge.version <- merge(ratio, response.meta, by = "sample", all.x = TRUE)
rownames(ratio) <- ratio$sample
ratio <- ratio[merge.version$sample ,]
ratio <- ratio[, -1]
ratio[is.na(ratio)] <- 0

# normalization
scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
scale_ratio <- as.data.frame(scale_ratio)
scale_ratio <- t(scale_ratio)



###绘图####

estim.coad <- readRDS(file.path(PROJECT_ROOT, "泛癌S", "NMF", "基质免疫样本2", "nmf_result.rds"))
plot(estim.coad)


nmf.rank5 = readRDS(file.path(PROJECT_ROOT, "泛癌S", "NMF", "基质免疫样本2", "nmf.rank6.rds"))
plot(nmf.rank5)

index <- extractFeatures(nmf.rank5, "max")  # change the order of the index

new.index <- list()
new.index[[1]] <- index[[1]]
new.index[[2]] <- index[[2]]
new.index[[3]] <- index[[4]]
new.index[[4]] <- index[[3]]
new.index[[5]] <- index[[5]]
new.index[[6]] <- index[[6]]

sig.order <- unlist(new.index)

NMF.Exp.rank5 <- scale_ratio[sig.order,]
NMF.Exp.rank5 <- na.omit(NMF.Exp.rank5)

group <- predict(nmf.rank5)  # adjust the position of the module

new.group <- c()
for (each in group) {
  if (each %in% c("1")) { new.group <- c(new.group, "1") }
  if (each %in% c("2")) { new.group <- c(new.group, "2") }
  if (each %in% c("4")) { new.group <- c(new.group, "3") }
  if (each %in% c("3")) { new.group <- c(new.group, "4") }
  if (each %in% c("5")) { new.group <- c(new.group, "5") }
  if (each %in% c("6")) { new.group <- c(new.group, "6") }
}

new.group <- factor(new.group, levels = c("1", "2", "3", "4","5","6"))#, "5"

z_ratio <- scale(ratio) / 4
z_ratio <- as.data.frame(z_ratio)
z_ratio <- t(z_ratio)

plot_matrix <- z_ratio[sig.order,]
plot_matrix <- na.omit(plot_matrix)

info.matrix <- as.data.frame(t(NMF.Exp.rank5))
info.matrix$sampleID <- rownames(info.matrix)
info.matrix$group <- new.group

colnames(info.matrix)[colnames(info.matrix) == "sampleID"] <- "sample"
info.matrix <- merge(info.matrix, response.meta, by = "sample", all.x = TRUE)

gene.group <- c()
for (each in rownames(NMF.Exp.rank5)) {
  if (each %in% rownames(scale_ratio)[new.index[[1]]]) {
    gene.group <- c(gene.group, "module1")
  } else if (each %in% rownames(scale_ratio)[new.index[[2]]]) {
    gene.group <- c(gene.group, "module2")
  } else if (each %in% rownames(scale_ratio)[new.index[[3]]]) {
    gene.group <- c(gene.group, "module3")
  } else if (each %in% rownames(scale_ratio)[new.index[[4]]]) {
    gene.group <- c(gene.group, "module4")
  } else if (each %in% rownames(scale_ratio)[new.index[[5]]]) {
    gene.group <- c(gene.group, "module5")
  } else if (each %in% rownames(scale_ratio)[new.index[[6]]]) {
    gene.group <- c(gene.group, "module6")
  } 
}



info.matrix[is.na(info.matrix)] <- "unknown"




library(dplyr)
library(stringr)

info.matrix <- info.matrix %>%
  mutate(
    treatment_zzm = {
      x <- str_squish(str_to_lower(as.character(treatments)))  # 统一大小写/空格
      case_when(
        is.na(x) ~ "unknow",                                                   # 真NA → unknow
        x %in% c("unknown","unknow","nan","na","none","no treatment") ~ "Unknow",
        x %in% c("untreated","naive","healthy") ~ "Untreated",                     # 明确未治疗
        TRUE ~ "Treated"                                                           # 其它统归治疗
      )
    },
    treatment_zzm = factor(treatment_zzm, levels = c("Untreated","Treated","Unknow"))
  )

# 查看结果
t
# 创建 HeatmapAnnotation 对象
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
 
  module = factor(info.matrix$group.x, levels = c("1", "2", "3", "4","5","6")),
  col = list(
    group = c(
      "Autoimmune diseases" = '#F1BB72',
      "Inflammation"         = '#D6E7A3',
      "Normal"               = '#E5D2DD',
      "Normal adjacent tissue" =  '#57C3F3',
      "Precancerous condition" = '#E59CC4',
      "Tumor"                = '#53A85F',
      "Tumor_metastasis"     =  '#F3B1A0',
      "Tumor_blood"          =  '#476D87',
      "Blood_Tumor"          = '#E95C59' ),
    tissue_organ = c("adnexa of uterus" = '#E5D2DD',"Adrenal" ='#53A85F', "Bladder" ='#F1BB72', "Blood" = '#F3B1A0',
                     "Blood_vessel" = '#D6E7A3', "Bone" = '#57C3F3', "Bone_marrow" = '#476D87', "Brain" = '#E95C59',
                     "breast" = '#E59CC4', "Esophagus" = '#AB3282', "Intestine" = '#23452F', "Kidney" = '#BD956A',
                     "Liver" = '#8C549C', "Lung" = '#585658', "Lymph_node" ='#9FA3A8', "Nose" ='#E0D4CA',
                     "Oral_cavity" ='#5F3D69', "Other" = '#C5DEBA', "Ovary" = '#58A4C3', "prostate" ='#E4C755',
                     "Skin" = '#F7F398', "spleen" = '#AA9A59', "Stomach" = '#E63863', "Uterus" ='#E39A35'),
    
    
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
    
    module = c("1" = "#E37D6F", "2" ="#81C8D9", "3" = "#129982", "4" ="#7584A7","5" ="#F1B4A3","6" ="#C1E6F3" )#, "5" = "#F39B7FB2"
  ),
  simple_anno_size = unit(0.5, "cm")
)

# 创建右侧颜色分组
gene_group_colors <- c(
  "module1" = "#E37D6F", 
  "module2" = "#81C8D9", 
  "module3" = "#129982", 
  "module4" = "#7584A7", 
  "module5" = "#F1B4A3",
  "module6" = "#C1E6F3"
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
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)


hp <- Heatmap(
  plot_matrix,
  col = col_fun,
  name = "ratio",
  top_annotation = ha,
  right_annotation = ra,
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
  use_raster = FALSE,         # 见下一条说明
  rect_gp = grid::gpar(col = NA, lwd = 0),  # ← 关键：去掉每个格子的边框
  border = NA                                  # ← 取消整张热图外边框
)
draw(hp)

Cairo::CairoPDF(file.path(PROJECT_ROOT, "泛癌S", "NMF", "基质免疫样本2", "cell_heatmap6.pdf"), width = 20, height = 15)
draw(hp, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
# 使用 ggplot2 将热图保存为文件
library(ggplot2)
heatmap_grob <- grid::grid.grab()
heatmap_ggplot <- cowplot::ggdraw() + cowplot::draw_grob(heatmap_grob)

# 保存为 PNG 和 PDF 格式
ggsave(heatmap_ggplot, filename = file.path(PROJECT_ROOT, "泛癌S", "NMF", "基质免疫样本2", "cell_heatmap6.png"), width = 30, height = 20, dpi = 300, device = 'png', bg = '#FFFFFF')
ggsave(heatmap_ggplot, filename = file.path(PROJECT_ROOT, "泛癌S", "NMF", "基质免疫样本2", "cell_heatmap6.pdf"), width = 30, height = 20, dpi = 300, device = 'pdf', bg = '#FFFFFF')



####绘制每个module种的分组占比######
table(info.matrix$group.x, useNA = "ifany")
info.matrix <- info.matrix[!is.na(info.matrix$group.x), ]



library(dplyr)
library(ggplot2)
library(grid)

label_cutoff <- 10

proportion_data <- filtered %>%
  dplyr::select(group.x, group.y) %>%
  dplyr::filter(
    !is.na(group.x),
    !is.na(group.y),
    group.x != "",
    group.y != ""
  ) %>%
  dplyr::group_by(group.x, group.y) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(group.x) %>%
  dplyr::mutate(
    percentage = count / sum(count) * 100,
    label = ifelse(
      percentage >= label_cutoff,
      sprintf("%.1f%%", percentage),
      NA_character_
    )
  ) %>%
  dplyr::ungroup()

group_x_levels <- sort(
  unique(as.numeric(as.character(proportion_data$group.x))),
  na.last = NA
)

proportion_data$group.x <- factor(
  proportion_data$group.x,
  levels = rev(group_x_levels)
)

group_y_colors <- c(
  "Autoimmune diseases"    = "#F1BB72",
  "Blood_Tumor"            = "#E95C59",
  "Inflammation"           = "#D6E7A3",
  "Normal"                 = "#E5D2DD",
  "Normal adjacent tissue" = "#57C3F3",
  "Precancerous condition" = "#E59CC4",
  "Tumor"                  = "#53A85F",
  "Tumor_metastasis"       = "#F3B1A0",
  # "Tumor_blood"            = "#476D87"
)

group_y_order <- c(
  "Tumor_metastasis",
  "Tumor",
  # "Tumor_blood",
  "Precancerous condition",
  "Normal adjacent tissue",
  "Normal",
  "Inflammation",
  "Blood_Tumor",
  "Autoimmune diseases"
)

proportion_data$group.y <- factor(
  proportion_data$group.y,
  levels = group_y_order
)

p <- ggplot(
  proportion_data,
  aes(x = group.x, y = percentage, fill = group.y)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.86,
    color = "white",
    linewidth = 0.15
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5,
    fontface = "bold",
    color = "black",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = group_y_colors,
    breaks = group_y_order,
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = function(x) paste0("Module ", x)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    breaks = seq(0, 100, 20)
  ) +
  coord_flip() +
  labs(
    title = "Percentage of Group Categories within Each Module",
    x = "Modules",
    y = "Percentage (%)",
    fill = "Group"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = unit(0.18, "cm"),
    
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    
    legend.position = "right",
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.55, "cm"),
    legend.background = element_blank(),
    
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.margin = margin(8, 12, 8, 8)
  )

print(p)



####绘制每个module种的肿瘤占比######
library(dplyr)
library(ggplot2)

my48colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175', '#9C6B30', '#497E76', '#9E003A', '#D9D2E9', '#7A9D96', '#EFAFB7',
                '#4D5061', '#8A9B0F', '#FFD1DC', '#C6D166', '#8C564B', '#7B4173')

colnames(info.matrix)
table(info.matrix$cancer)
# Assuming info.matrix is already available from your provided code
# Calculate the proportion of group.y within each group.x
proportion_data <- info.matrix %>%
  group_by(group.x, cancer) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group.x) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# 已有的基础配色（8种）
group_y_colors <- c(
  "Autoimmune diseases" = "#E5D2DD",
  "Inflammation" = "#53A85F",
  "Normal" = "#F1BB72",
  "Normal adjacent tissue" = "#F3B1A0",
  "Precancerous condition" = "#D6E7A3",
  "Tumor" = "#57C3F3",
  "Blood_Tumor" = "#476D87",
  "Tumor_metastasis" = "#E95C59"
)

# 为每个疾病类型单独分配颜色（从 my48colors 中按顺序取）
group_y_colors <- c(
  group_y_colors,
  "ALC" = "#23452F",
  "atherosclerosis" = "#BD956A",
  "atrophic gastritis" = "#8C549C",
  "B-ALL" = "#585658",
  "BE" = "#9FA3A8",
  "BLCA" = "#E0D4CA",
  "BRCA" = "#5F3D69",
  "CD" = "#C5DEBA",
  "CESC" = "#58A4C3",
  "CF" = "#E4C755",
  "CLL" = "#F7F398",
  "COPD" = "#AA9A59",
  "COVID-19" = "#E63863",
  "CR" = "#E39A35",
  "CRC" = "#C1E6F3",
  "cSCC" = "#6778AE",
  "CTCL" = "#91D0BE",
  "ESCC" = "#B53E2B",
  "FL" = "#712820",
  "Gastritis" = "#DCC1DD",
  "GBM" = "#CCE0F5",
  "GIM" = "#CCC9E6",
  "HCC" = "#625D9E",
  "HLE" = "#68A180",
  "HNSCC" = "#3A6963",
  "ICC" = "#968175",
  "ILD" = "#9C6B30",
  "KIRC" = "#497E76",
  "KO" = "#9E003A",
  "LAM" = "#D9D2E9",
  "LBTC" = "#7A9D96",
  "LCC" = "#EFAFB7",
  "LCNEC" = "#4D5061",
  "light chain disease" = "#8A9B0F",
  "LUAD" = "#FFD1DC",
  "LUAD+LCNEC" = "#C6D166",
  "LUSC" = "#8C564B",
  "MM" = "#7B4173",
  "Normal" = "#F1BB72",  # 保留原始定义
  "NPC" = "#E59CC4",
  "NSCLC" = "#AB3282",
  "OS" = "#23452F",
  "OSCC" = "#BD956A",
  "OV" = "#8C549C",
  "PBC" = "#585658",
  "PC" = "#9FA3A8",
  "PDAC" = "#E0D4CA",
  "PF" = "#5F3D69",
  "PNA" = "#C5DEBA",
  "PRAD" = "#58A4C3",
  "PSC" = "#E4C755",
  "RA" = "#F7F398",
  "sAH" = "#AA9A59",
  "SCLC" = "#E63863",
  "SS" = "#E39A35",
  "STAD" = "#C1E6F3",
  "UC" = "#6778AE",
  "UCEC" = "#91D0BE"
)

# 检查总颜色数是否足够
if (length(unique(names(group_y_colors))) < length(all_diseases)) {
  warning("颜色数量不足，建议补充新颜色！")
}


# 你的备用 48 种颜色（排除已使用的颜色）
my48colors <- my48colors[!my48colors %in% group_y_colors]

# 需要配色的疾病类型（从 table(info.matrix$cancer) 中提取）
all_diseases <- c(
  "ALC", "atherosclerosis", "atrophic gastritis", "B-ALL", "BE", "BLCA", "BRCA",
  "CD", "CESC", "CF", "CLL", "COPD", "COVID-19", "CR", "CRC", "cSCC", "CTCL",
  "ESCC", "FL", "Gastritis", "GBM", "GIM", "HCC", "HLE", "HNSCC", "ICC", "ILD",
  "KIRC", "KO", "LAM", "LBTC", "LCC", "LCNEC", "light chain disease", "LUAD",
  "LUAD+LCNEC", "LUSC", "MM", "Normal", "NPC", "NSCLC", "OS", "OSCC", "OV",
  "PBC", "PC", "PDAC", "PF", "PNA", "PRAD", "PSC", "RA", "sAH", "SCLC", "SS",
  "STAD", "UC", "UCEC"
)


# Create the stacked bar plot with percentage labels
p <- ggplot(proportion_data, aes(x = group.x, y = percentage, fill = cancer)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_y_colors) +
  labs(
    title = "Percentage of group.y Categories within Each cancer Module",
    x = "Group X (Modules)",
    y = "Percentage (%)",
    fill = "cancer"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(p)














