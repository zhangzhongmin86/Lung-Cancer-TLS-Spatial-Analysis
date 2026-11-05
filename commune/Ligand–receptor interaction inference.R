#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(future)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

arg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("counts", "genes", "metadata", "output-dir")
missing <- required[!required %in% names(arg)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(paste0("--", missing), collapse = ", "))
}

counts_file <- normalizePath(arg[["counts"]], mustWork = TRUE)
genes_file <- normalizePath(arg[["genes"]], mustWork = TRUE)
metadata_file <- normalizePath(arg[["metadata"]], mustWork = TRUE)
out_dir <- arg[["output-dir"]]
min_cells <- as.integer(arg[["min-cells"]] %||% "3")
p_cutoff <- as.numeric(arg[["p-cutoff"]] %||% "0.05")
min_interactions <- as.integer(arg[["min-interactions"]] %||% "500")
max_celltypes <- as.integer(arg[["max-celltypes"]] %||% "0")
seed <- as.integer(arg[["seed"]] %||% "1")
nboot <- as.integer(arg[["nboot"]] %||% "100")
groups_requested <- if (!is.null(arg[["groups"]])) {
  strsplit(arg[["groups"]], ",", fixed = TRUE)[[1]]
} else {
  NULL
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
future::plan("sequential")
set.seed(seed)

message("Reading sparse count matrix")
connection <- if (grepl("[.]gz$", counts_file)) gzfile(counts_file, "rt") else counts_file
counts <- as(Matrix::readMM(connection), "CsparseMatrix")
if (inherits(connection, "connection")) close(connection)
genes <- readLines(genes_file, warn = FALSE)
metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(counts) != length(genes)) {
  stop("Gene file length does not match count-matrix rows.")
}
if (ncol(counts) != nrow(metadata)) {
  stop("Metadata rows do not match count-matrix columns.")
}
required_columns <- c("sample", "celltype_3")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0) {
  stop("Metadata is missing columns: ", paste(missing_columns, collapse = ", "))
}
rownames(counts) <- make.unique(genes)
colnames(counts) <- rownames(metadata)

if (!is.null(groups_requested)) {
  if (!"group" %in% colnames(metadata)) stop("--groups requires a metadata group column")
  keep <- metadata$group %in% groups_requested
  metadata <- metadata[keep, , drop = FALSE]
  counts <- counts[, keep, drop = FALSE]
}

valid <- !is.na(metadata$sample) & !is.na(metadata$celltype_3) &
  metadata$celltype_3 != "" & metadata$celltype_3 != "nan"
metadata <- metadata[valid, , drop = FALSE]
counts <- counts[, valid, drop = FALSE]

library_size <- Matrix::colSums(counts)
scale_factor <- ifelse(library_size > 0, 10000 / library_size, 0)
normalized <- counts %*% Matrix::Diagonal(x = scale_factor)
normalized@x <- log1p(normalized@x)
dimnames(normalized) <- dimnames(counts)

samples <- sort(unique(as.character(metadata$sample)))
db <- CellChatDB.human
long_results <- list()
qc <- list()

for (sample_id in samples) {
  message("Running CellChat: ", sample_id)
  sample_cells <- rownames(metadata)[metadata$sample == sample_id]
  sample_meta <- metadata[sample_cells, , drop = FALSE]
  group_counts <- table(sample_meta$celltype_3)
  keep_groups <- names(group_counts)[group_counts >= min_cells]
  if (max_celltypes > 0 && length(keep_groups) > max_celltypes) {
    keep_groups <- names(sort(group_counts[keep_groups], decreasing = TRUE))[
      seq_len(max_celltypes)
    ]
  }
  sample_cells <- rownames(sample_meta)[sample_meta$celltype_3 %in% keep_groups]

  if (length(keep_groups) < 2 || length(sample_cells) == 0) {
    qc[[sample_id]] <- data.frame(
      sample = sample_id,
      cells = length(sample_cells),
      celltypes = length(keep_groups),
      interactions = 0,
      retained = FALSE,
      reason = "fewer than two eligible cell types"
    )
    next
  }

  sample_matrix <- normalized[, sample_cells, drop = FALSE]
  sample_meta <- sample_meta[sample_cells, , drop = FALSE]
  sample_meta$samples <- factor(sample_id)
  cellchat <- createCellChat(
    object = sample_matrix,
    meta = sample_meta,
    group.by = "celltype_3"
  )
  cellchat@DB <- db
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    population.size = TRUE,
    nboot = nboot,
    seed.use = seed
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  interactions <- subsetCommunication(cellchat, thresh = p_cutoff)

  if (is.null(interactions)) interactions <- data.frame()
  retained <- nrow(interactions) >= min_interactions
  qc[[sample_id]] <- data.frame(
    sample = sample_id,
    cells = length(sample_cells),
    celltypes = length(keep_groups),
    interactions = nrow(interactions),
    retained = retained,
    reason = if (retained) "retained" else "below min_interactions"
  )
  if (!retained) next

  interactions <- interactions %>%
    arrange(source, target, ligand, receptor, pval) %>%
    distinct(source, target, ligand, receptor, .keep_all = TRUE) %>%
    mutate(
      sample = sample_id,
      interaction_key = paste(source, target, ligand, receptor, sep = "|||"),
      .before = 1
    )
  long_results[[sample_id]] <- interactions
  message("  retained interactions: ", nrow(interactions))
}

qc_table <- bind_rows(qc)
write.csv(qc_table, file.path(out_dir, "cellchat_sample_qc.csv"), row.names = FALSE)
df_long <- bind_rows(long_results)
if (nrow(df_long) == 0) {
  stop(
    "No sample passed CellChat filtering. For a demo, lower --min-interactions; ",
    "for the manuscript analysis the intended threshold is 500."
  )
}
write.csv(
  df_long,
  file.path(out_dir, "cellchat_interactions_long.csv"),
  row.names = FALSE,
  quote = TRUE
)

interaction_keys <- sort(unique(df_long$interaction_key))
retained_samples <- sort(unique(df_long$sample))
probability <- matrix(
  0,
  nrow = length(interaction_keys),
  ncol = length(retained_samples),
  dimnames = list(interaction_keys, retained_samples)
)
p_value <- matrix(
  NA_real_,
  nrow = length(interaction_keys),
  ncol = length(retained_samples),
  dimnames = list(interaction_keys, retained_samples)
)

for (sample_id in retained_samples) {
  sample_data <- df_long[df_long$sample == sample_id, , drop = FALSE]
  probability[sample_data$interaction_key, sample_id] <- sample_data$prob
  p_value[sample_data$interaction_key, sample_id] <- sample_data$pval
}

write.csv(
  probability,
  file.path(out_dir, "cellchat_probability_interaction_by_sample.csv"),
  quote = TRUE
)
write.csv(
  p_value,
  file.path(out_dir, "cellchat_pvalue_interaction_by_sample.csv"),
  quote = TRUE,
  na = ""
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message(
  "Saved CellChat matrix: ", nrow(probability), " interactions x ",
  ncol(probability), " samples"
)
