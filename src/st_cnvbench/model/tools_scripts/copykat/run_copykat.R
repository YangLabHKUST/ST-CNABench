#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript run_copykat.R <input_10x_dir> <output_dir> <normal_cells_file> <genome> <win_size> <n_cores> <n_clust> <sample_name>")
}

# Arguments
input_dir   <- args[1] 
output_dir  <- args[2]
norm_file   <- args[3] 
genome_ver  <- args[4]
win_size    <- as.numeric(args[5])
n_cores     <- as.integer(args[6])
n_clust     <- as.integer(args[7])
sample_name <- args[8]

setwd(output_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(data.table)
})

message(paste0(">>> [CopyKAT] Step 1: Loading data..."))
# MTX
raw.data <- Read10X(data.dir = input_dir, gene.column = 1)
raw.data <- as.matrix(raw.data)

message(paste0("    Dimensions: ", nrow(raw.data), " genes x ", ncol(raw.data), " cells"))
if (norm_file != "None" && file.exists(norm_file)) {
  message(paste0("    Normal reference file: ", norm_file))
  norm.cell.names <- c()
  norm_df <- fread(norm_file, header = FALSE)
  norm.cell.names <- norm_df$V1
  norm.cell.names <- intersect(norm.cell.names, colnames(raw.data))
  message(paste0(">>> [CopyKAT] Found ", length(norm.cell.names), " valid normal reference cells."))

  message(">>> [CopyKAT] Step 2: Running main algorithm...")
  res <- copykat(
    rawmat = raw.data,
    id.type = "E",
    ngene.chr = 5, # modify to contain more chr
    win.size = win_size,
    sam.name = sample_name,
    distance = "euclidean",
    norm.cell.names = norm.cell.names,
    output.seg = FALSE,  
    plot.genes = TRUE, 
    genome = genome_ver,
    n.cores = n_cores
  )
} else {
  message(">>> [CopyKAT] No normal reference cells provided. Running in blind mode.")
  message(">>> [CopyKAT] Step 2: Running main algorithm...")
  res <- copykat(
    rawmat = raw.data,
    id.type = "E", # Use ensembl IDs
    ngene.chr = 5, # modify to contain more chr
    win.size = win_size,
    sam.name = sample_name,
    distance = "euclidean",
    #norm.cell.names = NULL,
    output.seg = FALSE,  
    plot.genes = TRUE, 
    genome = genome_ver,
    n.cores = n_cores
  )
}

# Extract Subclone predictions (<sample_name>_copykat_clustering_results.rds)
subclone_res <- readRDS(paste0(sample_name, "_copykat_clustering_results.rds"))
subclone <- cutree(subclone_res, k = n_clust)
# replace AGCT.1 to AGCT-1
names(subclone) <- gsub("\\.", "-", names(subclone))
df <- data.frame(
  Barcodes = names(subclone),
  Label_preds = paste0("copykat_", subclone),
  stringsAsFactors = FALSE
)
write.table(df, file = file.path(output_dir, "copykat_subcluster_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
message(paste0(">>> [CopyKAT] Subcluster predictions saved to: ", output_dir, "/copykat_subcluster_results.txt"))
# Save
saveRDS(res, file = file.path(output_dir, "copykat_object.rds"))
message(">>> [CopyKAT] Done.")