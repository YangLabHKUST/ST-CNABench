#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript run_numbat.R <matrix_dir> <output_dir> <meta_file> <allele_file> <n_clones> <n_cores> <genome_ver>")
}

# Arguments
matrix_dir  <- args[1]
output_dir  <- args[2]
meta_file   <- args[3]
allele_file <- args[4]
n_clones    <- as.numeric(args[5])
n_cores     <- as.numeric(args[6])
genome_ver  <- args[7]

setwd(output_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(numbat)
})

message(">>> [Numbat] Step 1: Loading Data...")
count_mat <- Read10X(data.dir = matrix_dir, gene.column = 2)

df_allele <- fread(allele_file)

if (meta_file != "None"){
  meta <- fread(meta_file)
  if ("Barcode" %in% colnames(meta)) setnames(meta, "Barcode", "cell")
  if ("barcode" %in% colnames(meta)) setnames(meta, "barcode", "cell")
  if ("tumor_normal" %in% colnames(meta)) setnames(meta, "tumor_normal", "group_type")
  if ("celltype_coarse" %in% colnames(meta)) setnames(meta, "celltype_coarse", "group") 

  # intersect cells
  common <- intersect(colnames(count_mat), meta$cell)
  message(sprintf(">>> [Numbat] Common barcodes: %d (Matrix: %d, Meta: %d)", 
                  length(common), ncol(count_mat), nrow(meta)))

  if (length(common) == 0) stop("No common cells found between matrix and metadata!")

  count_mat <- count_mat[, common, drop = FALSE]
  meta <- meta[cell %in% common]

  message(">>> [Numbat] Step 2: Aggregating Reference Counts...")
  ref_cells <- meta[group_type == "normal" | group_type == "Normal", cell]
  message(paste0("    Reference cells found: ", length(ref_cells)))

  if (length(ref_cells) == 0) stop("No normal cells found in metadata! Check 'tumor_normal' column.")

  # if no "group" column, assign all normal
  if (!"group" %in% colnames(meta)) {
      meta$group <- "Normal"
  }
  ref_annot <- meta[cell %in% ref_cells, .(cell, group)]
  setkey(ref_annot, cell)
  ref_annot <- ref_annot[ref_cells]
  lambdas_ref <- aggregate_counts(
    count_mat[, ref_cells, drop = FALSE],
    ref_annot,
    normalized = TRUE,
    verbose = TRUE
  )

  target_cells <- setdiff(colnames(count_mat), ref_cells)
  count_mat_target <- count_mat[, target_cells, drop = FALSE]
  message(paste0("    Target cells for analysis: ", length(target_cells)))

} else {
  message(">>> [Numbat] Step 2: No metadata file provided. Running in no reference normal mode.")
  # Maybe this ref_hca is not suitable for spatial data because the results find only 6 normal cells
  lambdas_ref <- numbat::ref_hca
  count_mat_target <- count_mat
}

message(">>> [Numbat] Step 3: Running Main Algorithm...")
out <- run_numbat(
    count_mat = count_mat_target,
    lambdas_ref = lambdas_ref,
    df_allele = df_allele,
    genome = genome_ver,
    t = 1e-5,
    max_entropy=0.8, # Replace max_entropy=0.5
    ncores = n_cores,
    plot = FALSE,
    n_cut = (n_clones-1),
    out_dir = output_dir
)

message(">>> [Numbat] Done.")