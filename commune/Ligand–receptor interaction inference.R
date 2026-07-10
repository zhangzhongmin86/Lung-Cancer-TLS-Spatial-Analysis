# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)




library(CellChat)
library(patchwork)


packageVersion("CellChat")



pbmc = readRDS("/data/beifen/zhongmin/泛癌S/h5ad/189数据/rds文件/数据/adata_part_1.rds")
table(pbmc$cancer)
table(pbmc$celltype_3)
zzm = pbmc@meta.data


table(pbmc$sample)


# 4) 简单检查：每个样本细胞数、celltype_3 分布
print(table(pbmc$sample))
print(table(pbmc$celltype_3))



suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
})


# 0) 基本检查

stopifnot("sample" %in% colnames(pbmc@meta.data))
stopifnot("celltype_3" %in% colnames(pbmc@meta.data))

samples <- sort(unique(as.character(pbmc$sample)))
print(table(pbmc$sample))

CellChatDB.use <- CellChatDB.human

min_cells <- 3
pv_cutoff <- 0.05  # 这里先保留，你现在没用到也没关系

# 存每个sample的df.net（长表）与 prob/pval 向量
df_list <- setNames(vector("list", length(samples)), samples)
prob_vec_list <- setNames(vector("list", length(samples)), samples)
pval_vec_list <- setNames(vector("list", length(samples)), samples)

set.seed(1)


# 1) 每个 sample 单独跑 CellChat

for (smp in samples) {
  message("Running CellChat for sample: ", smp)
  
  obj <- subset(pbmc, subset = sample == smp)
  
  meta <- obj@meta.data
  meta$celltype_3 <- as.character(meta$celltype_3)
  
  # 去掉 celltype_3 缺失的细胞（避免 CellChat 报错/产生NA群）
  keep_cells <- rownames(meta)[!is.na(meta$celltype_3) & meta$celltype_3 != ""]
  if (length(keep_cells) == 0) {
    message("  [skip] no valid cells after filtering celltype_3")
    next
  }
  obj <- subset(obj, cells = keep_cells)
  
  # 过滤每个 sample 内过小的 celltype（你设为3）
  group_counts <- table(obj$celltype_3)
  keep_groups <- names(group_counts)[group_counts >= min_cells]
  if (length(keep_groups) < 2) {
    message("  [skip] <2 celltype groups with >=", min_cells, " cells")
    next
  }
  obj <- subset(obj, subset = celltype_3 %in% keep_groups)
  
  # 表达矩阵（log-normalized）
  data.use <- GetAssayData(obj, assay = "RNA", slot = "data")
  meta.use <- obj@meta.data
  
  cellchat <- createCellChat(object = data.use, meta = meta.use, group.by = "celltype_3")
  cellchat <- setIdent(cellchat, ident.use = "celltype_3")
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  
  df.net <- subsetCommunication(cellchat)
  
  if (is.null(df.net) || nrow(df.net) == 0) {
    message("  [skip] df.net empty")
    next
  }
  
  # 生成唯一键：source_target_ligand_receptor
  df.net$interaction_key <- paste(df.net$source, df.net$target, df.net$ligand, df.net$receptor, sep = "_")
  
  # 若同一key重复，保留 pval 最小的一条（保险）
  df.net <- df.net %>%
    arrange(interaction_key, pval) %>%
    distinct(interaction_key, .keep_all = TRUE)
  
  df_list[[smp]] <- df.net
  
  pv <- df.net$prob
  names(pv) <- df.net$interaction_key
  prob_vec_list[[smp]] <- pv
  
  qv <- df.net$pval
  names(qv) <- df.net$interaction_key
  pval_vec_list[[smp]] <- qv
  
  message("  done: interactions=", nrow(df.net), ", groups=", length(levels(cellchat@idents)))
}

# 去掉没跑出结果的 sample
prob_vec_list <- prob_vec_list[!vapply(prob_vec_list, is.null, logical(1))]
pval_vec_list <- pval_vec_list[!vapply(pval_vec_list, is.null, logical(1))]
df_list <- df_list[!vapply(df_list, is.null, logical(1))]

if (length(prob_vec_list) == 0) stop("No sample produced valid CellChat results.")


# 2) 拼成矩阵：行=interaction_key，列=sample

all_keys <- sort(unique(unlist(lapply(prob_vec_list, names))))
sample_ok <- names(prob_vec_list)

prob_mat <- matrix(
  NA_real_,
  nrow = length(all_keys),
  ncol = length(sample_ok),
  dimnames = list(all_keys, sample_ok)
)

pval_mat <- matrix(
  NA_real_,
  nrow = length(all_keys),
  ncol = length(sample_ok),
  dimnames = list(all_keys, sample_ok)
)

for (smp in sample_ok) {
  pv <- prob_vec_list[[smp]]
  prob_mat[names(pv), smp] <- as.numeric(pv)
  
  qv <- pval_vec_list[[smp]]
  pval_mat[names(qv), smp] <- as.numeric(qv)
}


# 3) 输出：矩阵 + 长表

out_dir <- "/data/beifen/zhongmin/泛癌S/h5ad/189数据/rds文件/数据/cellchat_by_sample"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_prob_csv <- file.path(out_dir, "cellchat_prob_matrix_interaction_by_sample.csv")
out_pval_csv <- file.path(out_dir, "cellchat_pval_matrix_interaction_by_sample.csv")
out_long_csv <- file.path(out_dir, "cellchat_dfnet_long_all_samples.csv")

write.csv(prob_mat, out_prob_csv, quote = FALSE)
write.csv(pval_mat, out_pval_csv, quote = FALSE)

df_long <- bind_rows(lapply(names(df_list), function(smp) {
  d <- df_list[[smp]]
  d$sample <- smp
  d
}))
write.csv(df_long, out_long_csv, row.names = FALSE, quote = FALSE)

message("Saved:")
message("  ", out_prob_csv)
message("  ", out_pval_csv)
message("  ", out_long_csv)
message("Matrix dim (prob): ", nrow(prob_mat), " x ", ncol(prob_mat),
        "  [rows=interaction_key, cols=sample]")


