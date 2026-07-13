#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(NMF)
  library(dplyr)
  library(tidyr)
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
required <- c("matrix", "interactions", "output-dir")
missing <- required[!required %in% names(arg)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(paste0("--", missing), collapse = ", "))
}

matrix_file <- normalizePath(arg[["matrix"]], mustWork = TRUE)
interactions_file <- normalizePath(arg[["interactions"]], mustWork = TRUE)
out_dir <- arg[["output-dir"]]
min_samples <- as.integer(arg[["min-samples"]] %||% "3")
top_features <- as.integer(arg[["top-features"]] %||% "10000")
rank <- as.integer(arg[["rank"]] %||% "14")
nrun <- as.integer(arg[["nrun"]] %||% "200")
seed <- as.integer(arg[["seed"]] %||% "100")
survey_nrun <- as.integer(arg[["survey-nrun"]] %||% "0")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
probability <- as.matrix(read.csv(
  matrix_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
))
storage.mode(probability) <- "double"
probability[!is.finite(probability)] <- 0

interaction_metadata <- read.csv(
  interactions_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  arrange(interaction_key, pval) %>%
  distinct(interaction_key, .keep_all = TRUE)

hla_keys <- interaction_metadata %>%
  filter(
    grepl("(^|_)HLA-[A-F]", ligand) |
      grepl("(^|_)HLA-[A-F]", receptor)
  ) %>%
  pull(interaction_key)
probability <- probability[!rownames(probability) %in% hla_keys, , drop = FALSE]
probability <- probability[rowSums(probability > 0) >= min_samples, , drop = FALSE]
input_samples <- colnames(probability)
probability <- probability[rowSums(probability) > 0, colSums(probability) > 0, drop = FALSE]
if (nrow(probability) < 2 || ncol(probability) < 2) {
  stop("Too few non-zero interactions or samples remain for NMF.")
}

feature_variance <- apply(probability, 1, var)
feature_variance[!is.finite(feature_variance)] <- 0
keep_n <- min(top_features, sum(feature_variance > 0))
if (keep_n < 2) stop("Fewer than two variable interaction features remain.")
keep <- order(feature_variance, decreasing = TRUE)[seq_len(keep_n)]
probability <- probability[keep, , drop = FALSE]

minmax <- function(x) {
  range_x <- range(x, finite = TRUE)
  if (diff(range_x) == 0) return(rep(0, length(x)))
  (x - range_x[[1]]) / diff(range_x)
}
scaled <- t(apply(probability, 1, minmax))
scaled <- scaled[rowSums(scaled) > 0, , drop = FALSE]

max_rank <- min(nrow(scaled), ncol(scaled)) - 1
if (rank > max_rank) {
  warning("Requested rank ", rank, " exceeds matrix limit; using ", max_rank)
  rank <- max_rank
}
if (rank < 2) stop("NMF rank must be at least 2 after matrix filtering.")

write.csv(probability, file.path(out_dir, "communication_matrix_filtered.csv"))
write.csv(scaled, file.path(out_dir, "communication_matrix_scaled.csv"))
write.csv(
  data.frame(
    sample = input_samples,
    retained = input_samples %in% colnames(scaled),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "nmf_sample_filter.csv"),
  row.names = FALSE
)

if (survey_nrun > 0 && max_rank >= 2) {
  survey_ranks <- 2:min(20, max_rank)
  set.seed(seed)
  survey <- nmf(scaled, survey_ranks, nrun = survey_nrun, method = "lee")
  saveRDS(survey, file.path(out_dir, "nmf_rank_survey.rds"))
  pdf(file.path(out_dir, "nmf_rank_survey.pdf"), width = 7, height = 5)
  plot(survey)
  dev.off()
}

set.seed(seed)
fit <- nmf(scaled, rank, nrun = nrun, method = "lee", seed = seed)
saveRDS(fit, file.path(out_dir, paste0("nmf_rank", rank, ".rds")))

w <- basis(fit)
h <- coef(fit)
feature_indices <- extractFeatures(fit, "max")
module_rows <- bind_rows(lapply(seq_along(feature_indices), function(module_id) {
  index <- feature_indices[[module_id]]
  if (is.character(index)) index <- match(index, rownames(w))
  index <- as.integer(index[!is.na(index)])
  if (length(index) == 0) return(NULL)
  data.frame(
    interaction_key = rownames(w)[index],
    module = module_id,
    weight = w[index, module_id],
    stringsAsFactors = FALSE
  )
})) %>%
  left_join(interaction_metadata, by = "interaction_key") %>%
  arrange(module, desc(weight))

sample_assignments <- data.frame(
  sample = colnames(h),
  module = apply(h, 2, which.max),
  coefficient = apply(h, 2, max),
  stringsAsFactors = FALSE
)

write.csv(
  module_rows,
  file.path(out_dir, "nmf_module_interactions.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  sample_assignments,
  file.path(out_dir, "nmf_sample_assignments.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(w, file.path(out_dir, "nmf_basis.csv"))
write.csv(h, file.path(out_dir, "nmf_coefficients.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message(
  "Saved rank-", rank, " NMF with ", nrow(module_rows),
  " module-specific interactions and ", ncol(h), " samples"
)
