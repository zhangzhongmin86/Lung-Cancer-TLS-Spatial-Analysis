#!/usr/bin/env Rscript

required <- c(
  CellChat = "2.2.0",
  Matrix = "1.7.3",
  NMF = "0.28",
  dplyr = "1.1.4",
  tidyr = "1.3.1",
  tibble = "3.2.1",
  future = "1.34.0"
)

failed <- FALSE
for (package in names(required)) {
  if (!requireNamespace(package, quietly = TRUE)) {
    cat(package, "MISSING\n")
    failed <- TRUE
  } else {
    cat(package, as.character(packageVersion(package)), "\n")
  }
}
if (failed) quit(status = 1)
