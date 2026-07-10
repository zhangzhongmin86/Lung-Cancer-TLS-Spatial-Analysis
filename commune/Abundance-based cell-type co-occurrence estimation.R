# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)


library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggpubr)
library(stringr)
library(tibble)

# ---- 读数据
obs_file <- "/data/beifen/zhongmin/泛癌S/h5ad/189数据/combined_adata_innerobs.csv"
metadata <- read.csv(obs_file, row.names = 1, check.names = TRUE)

# ---- 免疫细胞定义（用 celltype_1 判断）
immune_ct1 <- c("T cells","Monocyte_macrophage","B cells","NK cells",
                "Plasma cells","DC","Mast cells","Neutrophils","pDC")

# ---- 只保留“既有免疫、又有非免疫”的样本
sample_ok <- metadata %>%
  group_by(sample) %>%
  summarise(
    has_immune   = any(celltype_1 %in% immune_ct1),
    has_nonimmune= any(!(celltype_1 %in% immune_ct1)),
    .groups = "drop"
  ) %>%
  filter(has_immune & has_nonimmune) %>%
  pull(sample)

dat <- metadata %>% filter(sample %in% sample_ok)

# ---- 输出目录
out_dir <- "/data/beifen/zhongmin/泛癌S/h5ad/癌和癌旁比较细胞类型/计算细胞之间的相关性"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 相关性分级（用于 Correlation_type 字段）
corr_label <- function(r){
  case_when(
    r >= 0.70 ~ "Strong positive",
    r >= 0.50 ~ "Moderate positive",
    r >= 0.30 ~ "Weak positive",
    r <= -0.70 ~ "Strong negative",
    r <= -0.50 ~ "Moderate negative",
    r <= -0.30 ~ "Weak negative",
    TRUE ~ "None"
  )
}

# ---- 按组循环计算
all_groups <- sort(unique(dat$group))

for (g in sort(unique(dat$group))) {
  df_g <- dat %>%
    dplyr::filter(group == g) %>%
    dplyr::mutate(
      sample    = as.character(sample),
      celltype_3 = as.character(celltype_3)
    )
  
  # 按样本计总细胞数（分母）
  denom <- dplyr::count(df_g, sample, name = "total_cells")
  
  # 每样本 x celltype_3 计数
  cnt <- dplyr::count(df_g, sample, celltype_3, name = "n")
  
  # 样本内比例 = 该 celltype_3 数量 / 该样本总细胞数
  prop_wide <- cnt %>%
    dplyr::left_join(denom, by = "sample") %>%
    dplyr::mutate(prop = n / total_cells) %>%
    dplyr::select(sample, celltype_3, prop) %>%
    tidyr::pivot_wider(names_from = celltype_3, values_from = prop, values_fill = 0) %>%
    dplyr::arrange(sample)
  
  if (ncol(prop_wide) <= 2) {
    message(sprintf("[skip] %s: 有效细胞类型列太少，未生成结果。", g))
    next
  }
  
  mat <- as.matrix(prop_wide[, -1, drop = FALSE])
  
  keep <- apply(mat, 2, sd, na.rm = TRUE) > 0
  mat  <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 2) {
    message(sprintf("[skip] %s: 方差>0 的细胞类型 < 2，未生成结果。", g))
    next
  }
  
  ct_names <- colnames(mat)
  pairs <- combn(ct_names, 2, simplify = FALSE)
  
  res_list <- lapply(pairs, function(p) {
    x <- mat[, p[1]]; y <- mat[, p[2]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
    tibble::tibble(
      Cell_Type_1 = p[1],
      Cell_Type_2 = p[2],
      Correlation = as.numeric(ct$estimate),
      P_value     = as.numeric(ct$p.value),
      N           = sum(ok)
    )
  })
  
  res <- dplyr::bind_rows(res_list)
  if (nrow(res) == 0) {
    message(sprintf("[skip] %s: 无法计算有效相关。", g))
    next
  }
  
  res <- res %>%
    dplyr::mutate(
      Adjusted_p_value = p.adjust(P_value, method = "BH"),
      Correlation_type = case_when(
        Correlation >= 0.70 ~ "Strong positive",
        Correlation >= 0.50 ~ "Moderate positive",
        Correlation >= 0.30 ~ "Weak positive",
        Correlation <= -0.70 ~ "Strong negative",
        Correlation <= -0.50 ~ "Moderate negative",
        Correlation <= -0.30 ~ "Weak negative",
        TRUE ~ "None"
      )
    ) %>%
    dplyr::filter(Correlation > 0.3, P_value < 0.05) %>%
    dplyr::select(Cell_Type_1, Cell_Type_2, Correlation, P_value, Adjusted_p_value, Correlation_type) %>%
    dplyr::arrange(dplyr::desc(Correlation), P_value)
  
  out_file <- file.path(out_dir, paste0("specificity_correlation_", g, ".csv"))
  write.csv(res, out_file, row.names = FALSE)
  message("Saved: ", out_file)
}

