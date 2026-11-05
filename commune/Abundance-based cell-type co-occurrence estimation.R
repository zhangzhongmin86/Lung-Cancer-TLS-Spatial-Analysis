#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
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
required <- c("metadata", "output-dir")
missing <- required[!required %in% names(arg)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(paste0("--", missing), collapse = ", "))
}

metadata_file <- normalizePath(arg[["metadata"]], mustWork = TRUE)
out_dir <- arg[["output-dir"]]
corr_cutoff <- as.numeric(arg[["correlation-cutoff"]] %||% "0.3")
adjusted_p_cutoff <- as.numeric(arg[["adjusted-p-cutoff"]] %||% "0.05")
min_samples <- as.integer(arg[["min-samples"]] %||% "3")
groups_requested <- if (!is.null(arg[["groups"]])) {
  strsplit(arg[["groups"]], ",", fixed = TRUE)[[1]]
} else {
  NULL
}

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_columns <- c("sample", "group", "celltype_1", "celltype_3")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0) {
  stop("Metadata is missing columns: ", paste(missing_columns, collapse = ", "))
}

immune_celltype_1 <- c(
  "T cells", "Myeloid cells", "B cells", "NK cells", "NKT cells", "Plasma cells"
)

sample_ok <- metadata %>%
  filter(!is.na(sample), !is.na(celltype_3), celltype_3 != "", celltype_3 != "nan") %>%
  group_by(sample) %>%
  summarise(
    has_immune = any(celltype_1 %in% immune_celltype_1),
    has_nonimmune = any(!(celltype_1 %in% immune_celltype_1)),
    .groups = "drop"
  ) %>%
  filter(has_immune & has_nonimmune) %>%
  pull(sample)

dat <- metadata %>%
  filter(sample %in% sample_ok, !is.na(group), !is.na(celltype_3), celltype_3 != "")
if (!is.null(groups_requested)) {
  dat <- dat %>% filter(group %in% groups_requested)
}
if (nrow(dat) == 0) {
  stop("No cells remain after immune/non-immune and group filtering.")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(
  out_dir,
  pattern = "^specificity_correlation_.*[.]csv$",
  full.names = TRUE
))
all_results <- list()

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

for (group_name in sort(unique(dat$group))) {
  group_data <- dat %>%
    filter(group == group_name) %>%
    transmute(sample = as.character(sample), celltype_3 = as.character(celltype_3))

  sample_count <- n_distinct(group_data$sample)
  if (sample_count < min_samples) {
    message("[skip] ", group_name, ": only ", sample_count, " samples")
    next
  }

  denominator <- count(group_data, sample, name = "total_cells")
  proportion_wide <- group_data %>%
    count(sample, celltype_3, name = "n") %>%
    left_join(denominator, by = "sample") %>%
    mutate(proportion = n / total_cells) %>%
    select(sample, celltype_3, proportion) %>%
    pivot_wider(names_from = celltype_3, values_from = proportion, values_fill = 0) %>%
    arrange(sample)

  abundance <- as.matrix(proportion_wide[, -1, drop = FALSE])
  abundance <- abundance[, apply(abundance, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(abundance) < 2) {
    message("[skip] ", group_name, ": fewer than two variable cell types")
    next
  }

  pairs <- combn(colnames(abundance), 2, simplify = FALSE)
  results <- bind_rows(lapply(pairs, function(pair) {
    x <- abundance[, pair[[1]]]
    y <- abundance[, pair[[2]]]
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < min_samples) return(NULL)
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = "pearson"))
    tibble(
      Cell_Type_1 = pair[[1]],
      Cell_Type_2 = pair[[2]],
      Correlation = unname(test$estimate),
      P_value = test$p.value,
      N = sum(keep)
    )
  }))

  if (nrow(results) == 0) next
  results <- results %>%
    mutate(
      Adjusted_p_value = p.adjust(P_value, method = "BH"),
      Correlation_type = case_when(
        Correlation >= 0.70 ~ "Strong positive",
        Correlation >= 0.50 ~ "Moderate positive",
        Correlation >= 0.30 ~ "Weak positive",
        Correlation <= -0.70 ~ "Strong negative",
        Correlation <= -0.50 ~ "Moderate negative",
        Correlation <= -0.30 ~ "Weak negative",
        TRUE ~ "None"
      ),
      group = group_name
    ) %>%
    filter(Correlation > corr_cutoff, Adjusted_p_value < adjusted_p_cutoff) %>%
    arrange(desc(Correlation), Adjusted_p_value)

  output_file <- file.path(
    out_dir,
    paste0("specificity_correlation_", safe_name(group_name), ".csv")
  )
  write.csv(results %>% select(-group), output_file, row.names = FALSE, quote = TRUE)
  all_results[[group_name]] <- results
  message("Saved ", nrow(results), " supported pairs: ", output_file)
}

combined <- bind_rows(all_results)
write.csv(combined, file.path(out_dir, "cooccurrence_all_groups.csv"), row.names = FALSE)
writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)
