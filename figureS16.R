# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)







library(AUCell)
library(ggplot2)

pbmc <- readRDS("/work/zhongmin/数据/泛癌S/转h5ad/上皮细胞/上皮细胞去除批次效应后_去除ZZM1_最小样本数50_与总的400万细胞取交集_选取对应癌和正常上皮.rds")
zzm <- pbmc@meta.data

table(pbmc$celltype_3)
table(pbmc$tissue_organ)

out_dir <- "/work/zhongmin/数据/泛癌S/转h5ad/上皮细胞/figureS16_cell_death_scores"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obs_file <- "/work/zhongmin/数据/泛癌S/转h5ad/上皮细胞/细胞死亡基因.csv"
gene <- read.csv(obs_file, check.names = FALSE)

fix_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & x != "" & x != "NA"]
  x <- toupper(x)
  alias <- c(
    NRF2 = "NFE2L2",
    PARP = "PARP1",
    AMPK = "PRKAA1",
    PI3KC3 = "PIK3C3",
    ROCK = "ROCK1",
    JNK = "MAPK8",
    GLNA = "GLUL",
    HSP90A = "HSP90AA1",
    TNFSF6 = "FASLG",
    TNFRSF6 = "FAS"
  )
  x[x %in% names(alias)] <- alias[x[x %in% names(alias)]]
  unique(x)
}

pathway_cols <- setdiff(colnames(gene), "No.")
raw_gene_sets <- lapply(pathway_cols, function(col) fix_gene_symbol(gene[[col]]))
names(raw_gene_sets) <- make.names(pathway_cols)

anti_death_genes <- list(
  Apoptosis = c("BCL2", "BCL2L1", "BCL2L2", "MCL1", "BCL2A1", "CFLAR", "BIRC2", "BIRC3", "BIRC5", "BIRC6", "XIAP", "AKT1", "RELA", "NFKB1", "TNFAIP3"),
  Pyroptosis = c("GPX4", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", "CHMP6", "CHMP7"),
  Ferroptosis = c("GPX4", "SLC7A11", "SLC3A2", "AIFM2", "FSP1", "NFE2L2", "GCLC", "GCLM", "FTH1", "FTL", "SLC40A1", "SQSTM1", "TNFAIP3"),
  Autophagy = c("MTOR", "AKT1", "BCL2", "BCL2L1", "RUBCN"),
  Necroptosis = c("CASP8", "CFLAR", "BIRC2", "BIRC3", "XIAP", "TNFAIP3", "FTH1", "FTL", "SQSTM1", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", "CHMP6", "CHMP7"),
  Cuproptosis = c("MTF1", "GLS", "CDKN2A", "ATP7A", "ATP7B"),
  Parthanatos = c("ARH3", "RNF146", "ADPRHL2", "OGG1"),
  Entotic.cell.death = character(0),
  Lysosome.dependent.cell.death = c("FTH1", "FTL", "SQSTM1"),
  Alkaliptosis = c("CA9"),
  Oxeiptosis = c("NFE2L2"),
  Disulfidptosis = character(0)
)
anti_death_genes <- lapply(anti_death_genes, fix_gene_symbol)

geneSets <- unlist(lapply(names(raw_gene_sets), function(pathway) {
  genes <- raw_gene_sets[[pathway]]
  anti <- intersect(genes, anti_death_genes[[pathway]])
  pro <- setdiff(genes, anti)
  sets <- list()
  sets[[paste0(pathway, "_pro")]] <- pro
  if (length(anti) > 0) sets[[paste0(pathway, "_anti")]] <- anti
  sets
}), recursive = FALSE)

exprMatrix <- pbmc@assays$RNA@data
rownames(exprMatrix) <- Features(pbmc)
colnames(exprMatrix) <- Cells(pbmc)

geneSets <- lapply(geneSets, function(x) intersect(x, rownames(exprMatrix)))
geneSets <- geneSets[lengths(geneSets) > 0]

gene_set_table <- do.call(rbind, lapply(names(geneSets), function(nm) {
  data.frame(
    gene_set = nm,
    pathway = sub("_(pro|anti)$", "", nm),
    direction = sub("^.*_(pro|anti)$", "\\1", nm),
    gene = geneSets[[nm]],
    stringsAsFactors = FALSE
  )
}))
write.csv(gene_set_table, file.path(out_dir, "cell_death_gene_sets_pro_anti_used.csv"), row.names = FALSE)

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = TRUE, nCores = 64)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)

auc_matrix <- getAUC(cells_AUC)
auc_df <- as.data.frame(t(auc_matrix))

sample_col <- "sample"
celltype_col <- "celltype_3"
organ_col <- "tissue_organ"

if (!all(c(sample_col, celltype_col, organ_col) %in% colnames(zzm))) {
  stop("pbmc@meta.data must contain sample, celltype_3 and tissue_organ.")
}

meta_use <- zzm[rownames(auc_df), c(sample_col, celltype_col, organ_col), drop = FALSE]
colnames(meta_use) <- c("sample", "celltype_3", "tissue_organ")
meta_use$sample <- as.character(meta_use$sample)
meta_use$celltype_3 <- as.character(meta_use$celltype_3)
meta_use$tissue_organ <- as.character(meta_use$tissue_organ)
meta_use$condition <- ifelse(grepl("Malignant cells", meta_use$celltype_3, ignore.case = TRUE), "Tumor", "Normal")
meta_use$cancer <- meta_use$tissue_organ
meta_use <- meta_use[!is.na(meta_use$sample) & !is.na(meta_use$cancer) & meta_use$cancer != "Other", , drop = FALSE]

cell_count <- as.data.frame.matrix(table(meta_use$cancer, meta_use$condition))
sample_count <- as.data.frame.matrix(xtabs(~ cancer + condition, unique(meta_use[, c("sample", "cancer", "condition")])))
write.csv(cell_count, file.path(out_dir, "cell_death_tissue_condition_cell_counts.csv"))
write.csv(sample_count, file.path(out_dir, "cell_death_tissue_condition_sample_counts.csv"))

score_df <- cbind(meta_use, auc_df[rownames(meta_use), , drop = FALSE])
score_cols <- names(geneSets)

sample_score <- aggregate(
  score_df[, score_cols, drop = FALSE],
  by = list(sample = score_df$sample, cancer = score_df$cancer, condition = score_df$condition),
  FUN = mean,
  na.rm = TRUE
)
write.csv(sample_score, file.path(out_dir, "cell_death_AUCell_sample_mean_scores.csv"), row.names = FALSE)

compare_one <- function(df, score_col) {
  tumor <- df[df$condition == "Tumor", score_col]
  normal <- df[df$condition == "Normal", score_col]
  if (length(tumor) == 0 || length(normal) == 0) return(NULL)
  p_value <- if (length(tumor) >= 2 && length(normal) >= 2) {
    wilcox.test(tumor, normal)$p.value
  } else {
    NA_real_
  }
  data.frame(
    gene_set = score_col,
    tumor_median = median(tumor, na.rm = TRUE),
    normal_median = median(normal, na.rm = TRUE),
    delta_median = median(tumor, na.rm = TRUE) - median(normal, na.rm = TRUE),
    n_tumor_sample = length(tumor),
    n_normal_sample = length(normal),
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

comparison_list <- lapply(split(sample_score, sample_score$cancer), function(df) {
  out <- do.call(rbind, lapply(score_cols, function(score_col) compare_one(df, score_col)))
  if (is.null(out)) return(NULL)
  out$cancer <- unique(df$cancer)[1]
  out
})
comparison_list <- Filter(Negate(is.null), comparison_list)

if (length(comparison_list) == 0) {
  stop("No tissue_organ has both malignant and normal epithelial samples after excluding Other.")
}

comparison <- do.call(rbind, comparison_list)
comparison$padj <- p.adjust(comparison$p_value, method = "BH")
comparison$pathway <- sub("_(pro|anti)$", "", comparison$gene_set)
comparison$direction <- sub("^.*_(pro|anti)$", "\\1", comparison$gene_set)
comparison$neg_log10_padj <- -log10(pmax(comparison$padj, .Machine$double.xmin))
write.csv(comparison, file.path(out_dir, "cell_death_AUCell_tumor_vs_normal_comparison.csv"), row.names = FALSE)
write.csv(comparison[!is.na(comparison$padj) & comparison$padj < 0.05, ],
          file.path(out_dir, "cell_death_AUCell_tumor_vs_normal_comparison_FDR0.05.csv"),
          row.names = FALSE)

plot_dot_heatmap <- function(direction_use) {
  selected_organs <- c("Lung", "Intestine")
  df <- comparison[
    comparison$direction == direction_use &
      comparison$cancer %in% selected_organs &
      !is.na(comparison$padj) &
      comparison$padj < 0.05,
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0) {
    message("No significant ", direction_use, " results with FDR < 0.05.")
    return(invisible(NULL))
  }
  cancer_order <- selected_organs[selected_organs %in% unique(df$cancer)]
  pathway_order <- names(sort(tapply(df$delta_median, df$pathway, median, na.rm = TRUE)))
  df$cancer <- factor(df$cancer, levels = cancer_order)
  df$pathway_label <- gsub("\\.", " ", df$pathway)
  pathway_label_order <- gsub("\\.", " ", pathway_order)
  df$pathway_label <- factor(df$pathway_label, levels = rev(pathway_label_order))
  df$change_direction <- ifelse(df$delta_median >= 0, "Higher in tumor", "Lower in tumor")
  df$change_direction <- factor(df$change_direction, levels = c("Higher in tumor", "Lower in tumor"))
  title_text <- ifelse(direction_use == "pro", "Pro-death signatures", "Anti-death signatures")
  p <- ggplot(df, aes(x = cancer, y = pathway_label)) +
    geom_point(
      aes(fill = change_direction, size = pmin(neg_log10_padj, 10)),
      shape = 21,
      color = "#2A2A2A",
      stroke = 0.22,
      alpha = 1,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "Higher in tumor" = "#B65E2E",
        "Lower in tumor" = "#4C9DCB"
      ),
      name = "Direction"
    ) +
    scale_size_continuous(name = "-log10(FDR)", range = c(1.6, 8.5), breaks = c(2, 4, 6, 8, 10)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title_text) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      panel.grid.major = element_line(color = "#E8E2DC", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(face = "bold", size = 11, color = "#222222", margin = margin(b = 6)),
      axis.text.y = element_text(size = 10.5, color = "#222222"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 12, 8, 8)
    )
  ggsave(
    file.path(out_dir, paste0("cell_death_AUCell_", direction_use, "_tumor_vs_normal_dot_heatmap.pdf")),
    p,
    width = 6.2,
    height = max(4.8, length(unique(df$pathway)) * 0.36 + 1.6),
    useDingbats = FALSE
  )
}

plot_dot_heatmap("pro")
plot_dot_heatmap("anti")




gene_set_table_pro <- gene_set_table %>%
  dplyr::filter(direction == "pro")

write.csv(
  gene_set_table_pro,
  file.path(out_dir, "cell_death_gene_sets_pro_only_used.csv"),
  row.names = FALSE
)



