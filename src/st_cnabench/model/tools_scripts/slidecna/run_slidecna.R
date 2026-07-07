#!/usr/bin/env Rscript

# Run SlideCNA on one standardized benchmark dataset.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SlideCNA)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_arg) == 0) {
  stop("Unable to determine run_slidecna.R path for sourcing helper utilities.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
source(file.path(script_dir, "slidecna_export_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 16) {
  stop(
    "Usage: run_slidecna.R <matrix_dir> <metadata_file> <positions_file> ",
    "<gene_annot_file> <output_dir> <plot_dir> <dataset_id> <chrom_ord_csv> ",
    "<roll_mean_window> <avg_bead_per_bin> <spatial> <pos> <pos_k> <ex_k> ",
    "<use_GO_terms> <max_k_silhouette>",
    call. = FALSE
  )
}

matrix_dir <- args[[1]]
metadata_path <- args[[2]]
positions_path <- args[[3]]
gene_annot_path <- args[[4]]
output_directory <- args[[5]]
plot_directory <- args[[6]]
dataset_name <- args[[7]]
chrom_ord <- strsplit(args[[8]], ",", fixed = TRUE)[[1]]
roll_mean_window <- as.integer(args[[9]])
avg_bead_per_bin <- as.integer(args[[10]])
parse_cli_logical <- function(value) {
  # Parse Python/YAML-style booleans passed through the model wrapper.
  normalized <- tolower(as.character(value))
  if (normalized %in% c("true", "1", "yes")) {
    return(TRUE)
  }
  if (normalized %in% c("false", "0", "no")) {
    return(FALSE)
  }
  stop("Cannot parse logical argument: ", value, call. = FALSE)
}

spatial <- parse_cli_logical(args[[11]])
pos <- parse_cli_logical(args[[12]])
pos_k <- as.integer(args[[13]])
ex_k <- as.integer(args[[14]])
use_GO_terms <- parse_cli_logical(args[[15]])
max_k_silhouette <- as.integer(args[[16]])

matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
barcodes_path <- file.path(matrix_dir, "barcodes.tsv.gz")

require_file <- function(path) {
  # Stop immediately if a required input file is missing.
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

read_10x_counts <- function(matrix_path, features_path, barcodes_path, keep_genes) {
  # Read a 10x matrix and return a gene-by-barcode count data.frame.
  features <- fread(features_path, header = FALSE)
  barcodes <- fread(barcodes_path, header = FALSE)[[1]]
  counts_sparse <- readMM(matrix_path)

  if (nrow(features) != nrow(counts_sparse)) {
    stop("Feature count does not match matrix rows.", call. = FALSE)
  }
  if (length(barcodes) != ncol(counts_sparse)) {
    stop("Barcode count does not match matrix columns.", call. = FALSE)
  }

  gene_ids <- features[[1]]
  keep_idx <- which(gene_ids %in% keep_genes)
  if (length(keep_idx) == 0) {
    stop("No count genes matched gene_pos GENE values.", call. = FALSE)
  }

  counts_sparse <- counts_sparse[keep_idx, , drop = FALSE]
  rownames(counts_sparse) <- gene_ids[keep_idx]
  colnames(counts_sparse) <- barcodes
  as.data.frame(as.matrix(counts_sparse), check.names = FALSE)
}

make_beads_df <- function(metadata_path, positions_path, expected_barcodes) {
  # Merge tumor/normal labels with spatial coordinates for SlideCNA.
  metadata <- fread(metadata_path)
  positions <- fread(positions_path)

  required_metadata <- c("Barcode", "tumor_normal")
  required_positions <- c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  missing_positions <- setdiff(required_positions, names(positions))
  if (length(missing_metadata) > 0) {
    stop("Metadata missing columns: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  if (length(missing_positions) > 0) {
    stop("Spatial positions missing columns: ", paste(missing_positions, collapse = ", "), call. = FALSE)
  }

  beads_df <- merge(
    metadata[, .(bc = Barcode, cluster_type = tumor_normal)],
    positions[, .(bc = barcode, pos_x = pxl_col_in_fullres, pos_y = pxl_row_in_fullres)],
    by = "bc",
    all = FALSE
  )
  beads_df[, cluster_type := fifelse(tolower(cluster_type) == "normal", "Normal", "Malignant")]
  setkey(beads_df, bc)
  beads_df <- beads_df[expected_barcodes]

  if (anyNA(beads_df$bc) || nrow(beads_df) != length(expected_barcodes)) {
    stop("Metadata/spatial coordinates do not cover all matrix barcodes.", call. = FALSE)
  }
  if (!all(c("Normal", "Malignant") %in% unique(beads_df$cluster_type))) {
    stop("cluster_type must contain both Normal and Malignant labels.", call. = FALSE)
  }

  as.data.frame(beads_df)
}

make_gene_pos <- function(gene_annot_path, count_gene_ids, chrom_ord) {
  # Convert benchmark genome annotation into SlideCNA gene_pos format.
  annot <- fread(
    gene_annot_path,
    header = FALSE,
    col.names = c("chr", "start", "end", "gene_name", "gene_id", "arm", "chr_arm", "cytoband")
  )
  annot <- annot[gene_id %in% count_gene_ids]
  annot[, chr := as.character(chr)]
  annot[chr == "MT", chr := "M"]
  annot[!grepl("^chr", chr), chr := paste0("chr", chr)]
  annot <- annot[chr %in% chrom_ord]
  annot[, chr := factor(chr, levels = chrom_ord)]
  setorderv(annot, c("chr", "start", "end", "gene_id"))
  annot[, chr := as.character(chr)]
  annot <- annot[!duplicated(gene_id)]
  annot[, rel_gene_pos := seq_len(.N), by = chr]

  gene_pos <- annot[, .(GENE = gene_id, chr, start, end, rel_gene_pos)]
  if (nrow(gene_pos) < 1000) {
    stop("Too few genes matched genome annotation: ", nrow(gene_pos), call. = FALSE)
  }
  as.data.frame(gene_pos)
}

for (path in c(matrix_path, features_path, barcodes_path, metadata_path, positions_path, gene_annot_path)) {
  require_file(path)
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)

features <- fread(features_path, header = FALSE)
count_gene_ids <- features[[1]]
gene_pos <- make_gene_pos(gene_annot_path, count_gene_ids, chrom_ord)
counts <- read_10x_counts(matrix_path, features_path, barcodes_path, gene_pos$GENE)
beads_df <- make_beads_df(metadata_path, positions_path, colnames(counts))
gene_pos <- gene_pos[match(rownames(counts), gene_pos$GENE), ]

message("Prepared SlideCNA inputs for ", dataset_name, ":")
message("  counts: ", nrow(counts), " genes x ", ncol(counts), " beads")
message("  beads_df: ", nrow(beads_df), " beads")
message("  gene_pos: ", nrow(gene_pos), " genes")
message("  output_directory: ", output_directory)
message("  plot_directory: ", plot_directory)

run_slide_cna(
  counts = counts,
  beads_df = beads_df,
  gene_pos = gene_pos,
  output_directory = output_directory,
  plot_directory = plot_directory,
  spatial = spatial,
  roll_mean_window = roll_mean_window,
  avg_bead_per_bin = avg_bead_per_bin,
  pos = pos,
  pos_k = pos_k,
  ex_k = ex_k,
  max_k_silhouette = max_k_silhouette,
  use_GO_terms = use_GO_terms,
  chrom_ord = chrom_ord
)

export_result <- export_slidecna_malig_outputs(output_directory = output_directory, model_name = "SlideCNA")
message("Exported SlideCNA malignant benchmark outputs:")
message("  label_file: ", export_result$label_file)
message("  profile_file: ", export_result$profile_file)
message("  n_barcodes: ", export_result$n_barcodes)
message("  n_profile_rows: ", export_result$n_profile_rows)
