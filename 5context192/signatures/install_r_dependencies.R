#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript")
}
script_dir <- normalizePath(
  dirname(sub("^--file=", "", script_arg)),
  winslash = "/",
  mustWork = TRUE
)
local_library <- file.path(script_dir, ".Rlib")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = local_library)
}

required <- c("ICAMS", "mSigAct", "cosmicsig", "mSigTools")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  remotes::install_github(
    "steverozen/mSigAct",
    ref = "v3.0.1-branch",
    lib = local_library,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never"
  )
}

for (package in required) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Failed to install required R package: ", package)
  }
  cat(package, as.character(packageVersion(package)), "\n")
}
