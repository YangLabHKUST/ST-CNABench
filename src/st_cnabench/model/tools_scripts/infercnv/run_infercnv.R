#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript run_infercnv.R <counts_path> <annotation_file> <gene_order_file> <ref_norm> <output_dir> <cutoff> <n_threads> <k_obs_groups>")
}

# Arguments
counts_path     <- args[1]
annotation_file <- args[2]
gene_order_file <- args[3]
ref_norm       <- args[4]
output_dir      <- args[5]
cutoff          <- as.numeric(args[6])
n_threads       <- as.numeric(args[7])
k_obs_groups    <- as.numeric(args[8]) 


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
  library(data.table)
  library(Seurat)
})

message(paste0(">>> [InferCNV] Step 1: Loading Data..."))

raw_counts <- Read10X(data.dir = counts_path, gene.column = 1)

message(paste0(">>> [InferCNV] Matrix dimensions: ", nrow(raw_counts), " x ", ncol(raw_counts)))

message(">>> [InferCNV] Step 2: Creating InferCNV object...")
if (ref_norm == "False") {
  message(">>> [InferCNV] Running in no reference normal mode.")
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = raw_counts,
    annotations_file = annotation_file,
    gene_order_file = gene_order_file,
    delim = "\t",
    ref_group_names = NULL
  )
} else {
  message(">>> [InferCNV] Running with reference normal mode.")
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = raw_counts,
    annotations_file = annotation_file,
    gene_order_file = gene_order_file,
    delim = "\t",
    ref_group_names = c("Normal")
  )
}


message(">>> [InferCNV] Step 3: Running analysis...")

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff,
  out_dir = output_dir,
  cluster_by_groups = F,    
  k_obs_groups = k_obs_groups,
  HMM = T,
  HMM_type = "i6",
  denoise = TRUE,
  num_threads = n_threads,
  analysis_mode = 'subclusters'
)

save(infercnv_obj, file = file.path(output_dir, "infercnv_obj.rdata"))
message(">>> [InferCNV] Done.")