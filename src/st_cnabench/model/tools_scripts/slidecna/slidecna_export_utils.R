#!/usr/bin/env Rscript

# Export benchmark-facing SlideCNA malignant-only outputs from R-native objects.

suppressPackageStartupMessages({
  library(data.table)
})

export_slidecna_malig_outputs <- function(output_directory, model_name = "SlideCNA") {
  # Export barcode-level malignant labels and clone-level CNA profiles from SlideCNA RDS outputs.
  cna_rds <- file.path(output_directory, "cnv_data.rds")
  md_bin_rds <- file.path(output_directory, "md_bin.rds")

  if (!file.exists(cna_rds)) {
    stop("Required file not found: ", cna_rds, call. = FALSE)
  }
  if (!file.exists(md_bin_rds)) {
    stop("Required file not found: ", md_bin_rds, call. = FALSE)
  }

  slidecna_data <- readRDS(cna_rds)
  md_bin <- as.data.table(readRDS(md_bin_rds))

  required_cna <- c("malig", "hc_sub_malig")
  missing_cna <- setdiff(required_cna, names(slidecna_data))
  if (length(missing_cna) > 0) {
    stop("cnv_data.rds missing required entries: ", paste(missing_cna, collapse = ", "), call. = FALSE)
  }

  hc_sub_malig <- as.data.table(slidecna_data$hc_sub_malig)
  malig_long <- as.data.table(slidecna_data$malig)

  required_md <- c("bc", "bin_all")
  missing_md <- setdiff(required_md, names(md_bin))
  if (length(missing_md) > 0) {
    stop("md_bin.rds missing required columns: ", paste(missing_md, collapse = ", "), call. = FALSE)
  }

  required_hc <- c("clone", "variable")
  missing_hc <- setdiff(required_hc, names(hc_sub_malig))
  if (length(missing_hc) > 0) {
    stop("hc_sub_malig missing required columns: ", paste(missing_hc, collapse = ", "), call. = FALSE)
  }

  required_malig <- c("clone", "GENE", "chr", "start", "end", "value")
  missing_malig <- setdiff(required_malig, names(malig_long))
  if (length(missing_malig) > 0) {
    stop("malig entry missing required columns: ", paste(missing_malig, collapse = ", "), call. = FALSE)
  }

  hc_sub_malig[, clone := as.character(clone)]
  hc_sub_malig[, variable := as.integer(as.character(variable))]
  md_bin[, bin_all := as.integer(bin_all)]
  malig_long[, clone := as.character(clone)]

  unique_clones <- sort(unique(hc_sub_malig$clone))
  clone_map <- data.table(
    clone = unique_clones,
    Clone_ID = sprintf("%s_%d", model_name, seq_along(unique_clones))
  )

  label_dt <- merge(hc_sub_malig, clone_map, by = "clone", all.x = TRUE)
  label_dt <- merge(
    md_bin[, .(Barcodes = bc, variable = bin_all)],
    label_dt[, .(variable, Label_preds = Clone_ID)],
    by = "variable",
    all = FALSE
  )
  label_dt <- unique(label_dt[, .(Barcodes, Label_preds)])

  profile_dt <- merge(
    malig_long[, .(clone, GENE, chr, start, end, value)],
    clone_map,
    by = "clone",
    all.x = TRUE
  )
  profile_dt <- profile_dt[
    ,
    .(CN_linear = mean(as.numeric(value), na.rm = TRUE)),
    by = .(Clone_ID, Chromosome = chr, Start = as.integer(start), End = as.integer(end), Gene = GENE)
  ]
  profile_dt[, CN_Score := log2(pmax(CN_linear, 0.01))]
  profile_dt[, ID := paste(Chromosome, Start, End, Gene, sep = "_")]
  profile_dt[, LOH_Status := NA_real_]
  profile_dt <- profile_dt[, .(Clone_ID, Chromosome, Start, End, ID, CN_Score, LOH_Status)]

  label_out <- file.path(output_directory, "slidecna_malig_barcode_labels.tsv.gz")
  profile_out <- file.path(output_directory, "slidecna_malig_clone_profiles.tsv.gz")

  fwrite(label_dt, label_out, sep = "\t")
  fwrite(profile_dt, profile_out, sep = "\t")

  invisible(
    list(
      label_file = label_out,
      profile_file = profile_out,
      n_barcodes = nrow(label_dt),
      n_profile_rows = nrow(profile_dt)
    )
  )
}
