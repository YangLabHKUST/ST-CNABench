#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript run_scevan.R <matrix_dir> <output_dir> <meta_file> <sample_name> <n_cores>")
}

# Arguments
matrix_dir  <- args[1]
output_dir  <- args[2]
meta_file   <- args[3]
sample_name <- args[4]
n_cores     <- as.numeric(args[5])

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SCEVAN)
})

# Due to SCEVAN BUG, use output dir as working dir and output in ./output
# So the output will be SCEVAN/output/...
setwd(output_dir)

message(">>> [SCEVAN] Step 1: Loading Data...")

# MTX
mat_sparse <- readMM(file.path(matrix_dir, "matrix.mtx.gz"))
mat <- as.matrix(mat_sparse) 
features <- fread(file.path(matrix_dir, "features.tsv.gz"), header=FALSE)
barcodes <- fread(file.path(matrix_dir, "barcodes.tsv.gz"), header=FALSE)
rownames(mat) <- make.unique(features$V1) 
colnames(mat) <- barcodes$V1

# Annot
if (!is.null(meta_file) && meta_file != "None" && file.exists(meta_file)) {
  # Model with normal reference mode
meta <- fread(meta_file, header = TRUE)

# Intersect cells
common_cells <- intersect(colnames(mat), meta$Barcode)
if (length(common_cells) == 0) {
  stop("Error: No common barcodes found between Matrix and Metadata!")
}
message(sprintf(">>> [SCEVAN] Matched cells: %d", length(common_cells)))

mat_filtered  <- mat[, common_cells, drop = FALSE]
meta_filtered <- meta[Barcode %in% common_cells]
setkey(meta_filtered, Barcode)
meta_filtered <- meta_filtered[common_cells]
# normal or Normal
normal_barcodes <- meta_filtered[tumor_normal == "normal" | tumor_normal == "Normal", Barcode]
tumor_barcodes  <- meta_filtered[tumor_normal == "tumor" | tumor_normal == "Tumor",  Barcode]

message(sprintf(">>> [SCEVAN] Normal cells: %d", length(normal_barcodes)))
message(sprintf(">>> [SCEVAN] Tumor cells:  %d", length(tumor_barcodes)))

if (length(normal_barcodes) == 0) stop("No normal cells found (tumor_normal == 'normal'). SCEVAN requires normal reference.")

message(">>> [SCEVAN] Step 2: Running pipelineCNA...")

res_meta <- pipelineCNA(
  count_mtx = mat_filtered,
  sample    = sample_name,
  organism  = "human",   
  par_cores = n_cores,
  norm_cell = normal_barcodes,     
  FIXED_NORMAL_CELLS = TRUE,        
  SUBCLONES = TRUE,
  plotTree  = FALSE
)
} else {
  # Model without normal reference mode
  message(">>> [SCEVAN] Step 2: Running pipelineCNA without normal reference...")

  res_meta <- pipelineCNA(
    count_mtx = mat,
    sample    = sample_name,
    organism  = "human",   
    par_cores = n_cores,      
    SUBCLONES = TRUE,
    plotTree  = FALSE
  )
}
out_file <- file.path(output_dir, paste0(sample_name, "_SCEVAN_res_meta.tsv"))
data.table::fwrite(
  data.table::as.data.table(res_meta, keep.rownames = "Barcode"),
  out_file,
  sep = "\t"
)

message(">>> [SCEVAN] Done. Results saved to: ", out_file)