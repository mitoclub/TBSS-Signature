#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript")
}
script_dir <- normalizePath(
  dirname(sub("^--file=", "", script_arg)),
  winslash = "/",
  mustWork = TRUE
)
local_library <- file.path(script_dir, ".Rlib")
if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}

required <- c("ICAMS", "mSigAct", "cosmicsig")
missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Run install_r_dependencies.R first."
  )
}

branches <- c("low_Ts", "high_Ts", "high_minus_low_Ts")
input_dir <- file.path(script_dir, "input")
results_dir <- file.path(script_dir, "results")
work_dir <- file.path(script_dir, "work", "msigact")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

read_branch <- function(branch) {
  path <- file.path(input_dir, paste0(branch, "_samples.txt"))
  if (!file.exists(path)) {
    stop("Missing prepared input: ", path)
  }
  table <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(table) != 96L) {
    stop("Expected 96 mutation contexts in ", path)
  }
  names(table)[-1L] <- paste(branch, names(table)[-1L], sep = "__")
  table
}

branch_tables <- lapply(branches, read_branch)
mutation_types <- branch_tables[[1L]][[1L]]
for (table in branch_tables[-1L]) {
  if (!identical(table[[1L]], mutation_types)) {
    stop("Prepared branch matrices use different SBS96 row orders")
  }
}
combined <- branch_tables[[1L]]
for (table in branch_tables[-1L]) {
  combined <- cbind(combined, table[-1L])
}
combined_path <- file.path(input_dir, "samples_mSigAct.txt")
write.table(
  combined,
  combined_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
input_catalog <- ICAMS::ReadCatalog(file = combined_path)

sigs <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96
if (is.null(dim(sigs)) || nrow(sigs) != 96L) {
  stop("cosmicsig COSMIC v3.3 GRCh37 SBS96 database is unavailable")
}

signatures_to_exclude <- unique(c(
  "SBS32", "SBS11", "SBS25", "SBS31", "SBS35", "SBS86", "SBS87",
  "SBS90", "SBS22", "SBS88", "SBS27", "SBS43", "SBS45", "SBS46",
  "SBS47", "SBS48", "SBS49", "SBS50", "SBS51", "SBS52", "SBS53",
  "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60",
  "SBS95", "SBS9", "SBS84", "SBS85", "SBS4", "SBS29", "SBS92",
  "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"
))
relevant_sigs <- sigs[, !(colnames(sigs) %in% signatures_to_exclude), drop = FALSE]

run_assignment <- function(mode, signatures, priors) {
  output_dir <- file.path(work_dir, mode)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  priors <- priors[colnames(signatures)]
  if (anyNA(priors) || any(priors <= 0) || any(priors > 1)) {
    stop("Invalid signature-presence priors for mode ", mode)
  }
  cat(
    "mSigAct:", mode, "with", ncol(signatures), "candidate signatures and",
    ncol(input_catalog), "spectra\n"
  )
  mSigAct::MAPAssignActivity(
    spectra = input_catalog,
    sigs = signatures,
    sigs.presence.prop = priors,
    output.dir = output_dir,
    max.level = ncol(signatures) - 1L,
    max.subsets = 100000L,
    p.thresh = 0.05 / ncol(signatures),
    num.parallel.samples = 1L,
    mc.cores.per.sample = 1L,
    seed = 20260901L
  )
  invisible(output_dir)
}

read_exposures <- function(mode) {
  files <- list.files(
    file.path(work_dir, mode),
    pattern = "[.]exposure[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("No exposure files found for mSigAct mode ", mode)
  }
  rows <- lapply(files, function(path) {
    table <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(table) < 2L) {
      stop("Malformed exposure file: ", path)
    }
    data.frame(
      Mode = mode,
      Sample = names(table)[2L],
      Signature = as.character(table[[1L]]),
      Activity = as.numeric(table[[2L]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

trailing_args <- commandArgs(trailingOnly = TRUE)
aggregate_only <- "--aggregate-only" %in% trailing_args
requested_arg <- grep("^--modes=", trailing_args, value = TRUE)
requested_modes <- if (length(requested_arg) == 0L) {
  c("uniform", "sigprofiler_priors", "uniform_top6")
} else {
  strsplit(sub("^--modes=", "", requested_arg[[1L]]), ",", fixed = TRUE)[[1L]]
}
allowed_modes <- c("uniform", "sigprofiler_priors", "uniform_top6")
if (!all(requested_modes %in% allowed_modes)) {
  stop("Unknown --modes value. Allowed: ", paste(allowed_modes, collapse = ","))
}

if (!aggregate_only && (
  "uniform" %in% requested_modes || "uniform_top6" %in% requested_modes
)) {
  uniform_priors <- setNames(rep(1, ncol(relevant_sigs)), colnames(relevant_sigs))
  run_assignment("uniform", relevant_sigs, uniform_priors)
}

if (!aggregate_only && "sigprofiler_priors" %in% requested_modes) {
  priors_path <- file.path(results_dir, "sigprofiler_priors.csv")
  if (!file.exists(priors_path)) {
    stop("Missing SigProfiler priors: ", priors_path)
  }
  prior_table <- read.csv(priors_path, stringsAsFactors = FALSE)
  keep <- intersect(prior_table$Signature, colnames(relevant_sigs))
  if (length(keep) < 2L) {
    stop("Fewer than two SigProfiler prior signatures are available to mSigAct")
  }
  prior_signatures <- relevant_sigs[, colnames(relevant_sigs) %in% keep, drop = FALSE]
  prior_vector <- setNames(prior_table$Prior, prior_table$Signature)
  run_assignment("sigprofiler_priors", prior_signatures, prior_vector)
}

if (!aggregate_only && "uniform_top6" %in% requested_modes) {
  uniform_activities <- read_exposures("uniform")
  totals <- aggregate(Activity ~ Signature, uniform_activities, sum)
  totals <- totals[order(totals$Activity, decreasing = TRUE), ]
  totals <- totals[totals$Activity > 0, ]
  top_n <- min(6L, nrow(totals))
  if (top_n < 2L) {
    stop("Uniform mSigAct run identified fewer than two non-zero signatures")
  }
  empirical <- head(totals, top_n)
  empirical$Prior <- empirical$Activity / sum(empirical$Activity)
  write.csv(
    empirical,
    file.path(results_dir, "msigact_empirical_priors.csv"),
    row.names = FALSE
  )
  top_signatures <- relevant_sigs[
    , colnames(relevant_sigs) %in% empirical$Signature, drop = FALSE
  ]
  empirical_priors <- setNames(empirical$Prior, empirical$Signature)
  run_assignment("uniform_top6", top_signatures, empirical_priors)
}

completed_modes <- allowed_modes[
  vapply(
    allowed_modes,
    function(mode) length(list.files(
      file.path(work_dir, mode),
      pattern = "[.]exposure[.]csv$",
      recursive = TRUE
    )) > 0L,
    logical(1)
  )
]
activities <- do.call(rbind, lapply(completed_modes, read_exposures))
sample_parts <- strsplit(activities$Sample, "__", fixed = TRUE)
if (!all(lengths(sample_parts) == 3L)) {
  stop("Unexpected mSigAct sample names in exposure outputs")
}
activities$Branch <- vapply(sample_parts, `[[`, character(1), 1L)
activities$Gene <- vapply(sample_parts, `[[`, character(1), 2L)
activities$SpectrumType <- vapply(sample_parts, `[[`, character(1), 3L)
group_total <- ave(
  activities$Activity,
  activities$Mode,
  activities$Sample,
  FUN = sum
)
activities$ActivityFraction <- ifelse(group_total > 0, activities$Activity / group_total, 0)
write.csv(
  activities,
  file.path(results_dir, "msigact_activities.csv"),
  row.names = FALSE
)

distance_rows <- list()
for (mode in completed_modes) {
  files <- list.files(
    file.path(work_dir, mode),
    pattern = "[.]distances[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  for (path in files) {
    table <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    sample_name <- basename(dirname(path))
    sample_fields <- strsplit(sample_name, "__", fixed = TRUE)[[1L]]
    if (length(sample_fields) != 3L) {
      stop("Unexpected sample directory in mSigAct quality output: ", path)
    }
    table$Mode <- mode
    table$Sample <- sample_name
    table$Branch <- sample_fields[[1L]]
    table$Gene <- sample_fields[[2L]]
    table$SpectrumType <- sample_fields[[3L]]
    table$SourceFile <- paste(
      mode,
      basename(dirname(path)),
      basename(path),
      sep = "/"
    )
    distance_rows[[length(distance_rows) + 1L]] <- table
  }
}
if (length(distance_rows) > 0L) {
  distances <- do.call(rbind, distance_rows)
  write.csv(
    distances,
    file.path(results_dir, "msigact_quality.csv"),
    row.names = FALSE
  )
}

metadata <- data.frame(
  Tool = c("mSigAct", "ICAMS", "cosmicsig"),
  Version = c(
    as.character(packageVersion("mSigAct")),
    as.character(packageVersion("ICAMS")),
    as.character(packageVersion("cosmicsig"))
  ),
  COSMICVersion = "3.3",
  GenomeBuild = "GRCh37",
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  file.path(results_dir, "msigact_run_metadata.csv"),
  row.names = FALSE
)
cat(
  "Aggregated", length(unique(activities$Sample)), "spectra across modes:",
  paste(completed_modes, collapse = ", "), "\n"
)
