#!/usr/bin/env Rscript

# Export benchmark-facing SlideCNA malignant-only outputs for an existing result directory.

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_arg) == 0) {
  stop("Unable to determine export_slidecna_eval_outputs.R path for sourcing helper utilities.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
source(file.path(script_dir, "slidecna_export_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: export_slidecna_eval_outputs.R <slidecna_result_dir>", call. = FALSE)
}

result_dir <- normalizePath(args[[1]], mustWork = FALSE)
if (!dir.exists(result_dir)) {
  stop("Result directory not found: ", result_dir, call. = FALSE)
}

export_result <- export_slidecna_malig_outputs(output_directory = result_dir, model_name = "SlideCNA")
message("Exported SlideCNA malignant benchmark outputs:")
message("  label_file: ", export_result$label_file)
message("  profile_file: ", export_result$profile_file)
message("  n_barcodes: ", export_result$n_barcodes)
message("  n_profile_rows: ", export_result$n_profile_rows)
