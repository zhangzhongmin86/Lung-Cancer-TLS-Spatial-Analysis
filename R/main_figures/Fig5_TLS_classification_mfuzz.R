#!/usr/bin/env Rscript

# Cluster TLS reference-alignment profiles with Mfuzz.
#
# Usage:
#   Rscript Fig5_TLS_classification_mfuzz.R <pre_mfuzz.tsv> <output_dir>
#
# Optional environment variables:
#   TLS_MFUZZ_CLUSTERS  Number of fuzzy c-means clusters (default: 10)
#   TLS_MFUZZ_SEED      Random seed (default: 1)

options(stringsAsFactors = FALSE)

# Compatibility helper for packages built under R >= 4.4.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(left, right) if (is.null(left)) right else left
}

required_packages <- c("Biobase", "Mfuzz")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Biobase)
  library(Mfuzz)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Expected two arguments: <pre_mfuzz.tsv> <output_dir>. ",
    "Received ", length(args), "."
  )
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_clusters <- as.integer(Sys.getenv("TLS_MFUZZ_CLUSTERS", unset = "10"))
random_seed <- as.integer(Sys.getenv("TLS_MFUZZ_SEED", unset = "1"))
if (is.na(n_clusters) || n_clusters < 2) {
  stop("TLS_MFUZZ_CLUSTERS must be an integer >= 2.")
}
if (is.na(random_seed)) {
  stop("TLS_MFUZZ_SEED must be an integer.")
}

profiles <- read.delim(
  input_path,
  row.names = 1,
  check.names = FALSE
)
profiles <- as.matrix(profiles)
storage.mode(profiles) <- "double"

if (nrow(profiles) <= n_clusters) {
  stop("The number of TLS profiles must exceed the requested cluster count.")
}
if (ncol(profiles) < 2) {
  stop("Each TLS profile must contain at least two positions.")
}
if (any(!is.finite(profiles))) {
  stop("The Mfuzz input contains non-finite values.")
}

expression_set <- new("ExpressionSet", exprs = profiles)
fuzzifier <- mestimate(expression_set)

set.seed(random_seed)
fit <- mfuzz(
  expression_set,
  c = n_clusters,
  m = fuzzifier
)

profile_position <- seq_len(ncol(fit$centers)) - 1L
center_correlation <- apply(
  fit$centers,
  1,
  function(values) {
    suppressWarnings(
      cor(values, profile_position, method = "spearman", use = "complete.obs")
    )
  }
)
center_slope <- apply(
  fit$centers,
  1,
  function(values) unname(coef(lm(values ~ profile_position))[[2]])
)

# A non-negative center trend is compatible with the ordered reference axis;
# a negative center trend is classified as reference-divergent.
cluster_group <- ifelse(
  center_correlation >= 0,
  "Conforming",
  "Deviating"
)

cluster_definitions <- data.frame(
  cluster = seq_len(n_clusters),
  center_spearman_rho = center_correlation,
  center_slope = center_slope,
  Group = cluster_group,
  n_TLS = tabulate(fit$cluster, nbins = n_clusters)
)

trajectory_groups <- data.frame(
  cluster = as.integer(fit$cluster),
  max_membership = apply(fit$membership, 1, max),
  center_spearman_rho = center_correlation[fit$cluster],
  center_slope = center_slope[fit$cluster],
  Group = unname(cluster_group[fit$cluster]),
  row.names = rownames(profiles)
)

group_path <- file.path(output_dir, "trajectory_groups.tsv")
cluster_path <- file.path(output_dir, "mfuzz_cluster_definitions.tsv")
center_path <- file.path(output_dir, "mfuzz_cluster_centers.tsv")
plot_path <- file.path(output_dir, "mfuzz_clusters.pdf")
session_path <- file.path(output_dir, "mfuzz_session_info.txt")

write.table(
  trajectory_groups,
  group_path,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  cluster_definitions,
  cluster_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  fit$centers,
  center_path,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

pdf(plot_path, width = 10, height = 9)
mfuzz.plot(
  expression_set,
  fit,
  mfrow = c(2, 5),
  new.window = FALSE,
  time.labels = colnames(profiles)
)
dev.off()

capture.output(sessionInfo(), file = session_path)

message("TLS profiles: ", nrow(profiles))
message("Mfuzz clusters: ", n_clusters)
message("Random seed: ", random_seed)
message("Fuzzifier m: ", signif(fuzzifier, 6))
message("Group counts:")
print(table(trajectory_groups$Group))
message("Saved trajectory groups: ", group_path)
