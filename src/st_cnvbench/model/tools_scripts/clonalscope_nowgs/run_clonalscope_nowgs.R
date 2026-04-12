#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript run_clonalscope_nowgs.R <matrix_dir> <output_dir> <tumor_normal_file> <gene_coords_file> <aux_data_dir> <mincell>")
}

# Arguments
matrix_dir      <- args[1]
output_dir      <- args[2]
meta_file      <- args[3]
gene_coords_file<- args[4]
aux_dir         <- args[5]
mincell         <- as.numeric(args[6])


if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

suppressPackageStartupMessages({
  library(Clonalscope)
  library(Matrix)
  library(Seurat)
  library(data.table)
})

message(">>> [Clonalscope_NoWGS] Step 1: Loading data...")

# MTX
raw_counts <- Read10X(data.dir = matrix_dir)
matrix_obj <- raw_counts

feat_path <- file.path(matrix_dir, "features.tsv.gz")
feat <- read.delim(feat_path, header = FALSE, stringsAsFactors = FALSE)
features_df <- data.frame(
  gene_id = feat$V1,
  gene_name = feat$V2,
  feature_type = if (ncol(feat) >= 3) feat$V3 else "Gene Expression",
  stringsAsFactors = FALSE
)
barcodes <- data.frame(barcode = colnames(matrix_obj), stringsAsFactors = FALSE)
gene_symbols = feat[,2]
# Annot
if (meta_file == "None") {
  message(">>> [Clonalscope_NoWGS] No meta_file provided, will infer normal cells/spots from data.")
  initial_normal_spots = FindNormalReference(counts=matrix_obj, gene_symbols=gene_symbols,method="pca")
  celltype0 = barcodes; rownames(celltype0) = barcodes[,1]
  celltype0$type = "tumor"
  celltype0[initial_normal_spots,"type"] = "normal"
  message(sprintf(">>> [Clonalscope_NoWGS] Inferred %d normal cells/spots.", length(initial_normal_spots)))
} else {
  message(">>> [Clonalscope_NoWGS] Meta_file provided, loading cell type annotation...")
  # col1: barcode, col2: type (Tumor/Normal)
  celltype0 <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(celltype0) <- c("barcode", "type") 

  # Clonalscope require "type" column to be "tumor" or "normal"
  celltype0$type <- ifelse(grepl("Tumor|tumor", celltype0$type, ignore.case = TRUE), "tumor", "normal")

  # Intersect cells
  common_cells <- intersect(barcodes$barcode, celltype0$barcode)
  message(sprintf(">>> [Data Check] Kept %d common cells. Filtered out: %d from Matrix, %d from Annotation.", length(common_cells), nrow(barcodes) - length(common_cells), nrow(celltype0) - length(common_cells)))
  if (length(common_cells) == 0) stop("No common cells found between Matrix and Annotation!")

  matrix_obj <- matrix_obj[, common_cells]
  barcodes <- data.frame(barcode = colnames(matrix_obj), stringsAsFactors = FALSE)
  celltype0 <- celltype0[celltype0$barcode %in% common_cells, ]
}
# Aux Data
size_file <- file.path(aux_dir, "sizes.cellranger-GRCh38-1.0.0.txt")
cycle_file <- file.path(aux_dir, "cyclegenes.rds")
arm_file <- file.path(aux_dir, "cytoarm_table_hg38.txt")
bin_bed_file <- file.path(aux_dir, "hg38_200kb.windows.bed")

if (!file.exists(size_file)) stop("Missing sizes.cellranger-GRCh38-1.0.0.txt")
if (!file.exists(cycle_file)) stop("Missing cyclegenes.rds")

size <- read.table(size_file, stringsAsFactors = FALSE)
cyclegenes <- readRDS(cycle_file)
chrarm <- read.table(arm_file, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
bin_bed_200k <- read.table(bin_bed_file, header = FALSE, stringsAsFactors = FALSE)

# Bed file should contain columns: gene chr start end, like "chr1" rather than "1"
bed <- read.table(gene_coords_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# Clonalscope require bed as chr, start, end, geneID
bed <- bed[, c(2, 3, 4, 1)]
colnames(bed) <- c("chr", "start", "end", "geneID")

bed$chr <- paste0("chr", bed$chr)

message(">>> [Clonalscope_NoWGS] Step 2: Filtering features...")
Input_filtered <- FilterFeatures(
  mtx = matrix_obj,
  barcodes = barcodes,
  features = features_df,
  cyclegenes = cyclegenes
)

# Seg table
chrarm <- chrarm[order(as.numeric(chrarm[, 1]), as.numeric(chrarm[, 3])), ]
bin_bed <- chrarm[, -2]
seg_table_arm <- data.frame(
  chr = bin_bed[, 1],
  start = as.numeric(bin_bed[, 2]),
  end = as.numeric(bin_bed[, 3]),
  states = 1,
  length = as.numeric(bin_bed[, 3]) - as.numeric(bin_bed[, 2]),
  mean = 0, var = 0, Var1 = 1:nrow(bin_bed), Freq = 50000,
  chrr = paste0(bin_bed[, 1], ":", bin_bed[, 2]),
  stringsAsFactors = FALSE
)

#clustering_barcodes <- celltype0[celltype0$type == "tumor", 1]
print(table(celltype0$type))

message(">>> [Clonalscope_NoWGS] Step 3: Running Round 1 (Arm Level)...")
Cov_obj <- RunCovCluster(
  mtx = Input_filtered$mtx,
  barcodes = Input_filtered$barcodes,
  features = Input_filtered$features,
  bed = bed,
  celltype0 = celltype0,
  var_pt = 0.99,
  var_pt_ctrl = 0.99,
  include = "all",
  alpha_source = "all",
  ctrl_region = NULL,
  seg_table_filtered = seg_table_arm,
  size = size,
  dir_path = output_dir,
  breaks = 50,
  ngene_filter=10,
  prep_mode = "intersect",
  est_cap = 2,
  clust_mode = "all",
  mincell = mincell
)

saveRDS(Cov_obj, file.path(output_dir, "Cov_obj_chrarm.rds"))

message(">>> [Clonalscope_NoWGS] Step 4: Creating refined segmentation (NoWGS mode)...")

# Reassign clusters based on clustering cell * 0.01 # See: https://github.com/seasoncloud/Clonalscope/tree/main/samples/P5847/scRNA
result <- AssignCluster(Cov_obj$result_final$clustering2, mincell = mincell)
Zest   <- result$Zest

# Clonalscope cannot read ensembl GTF, thus use bed as pseudo-GTF
gtf <- data.frame(
  V1 = gsub("^chr", "", bed$chr),
  V4 = bed$start,
  V5 = bed$end,
  gene_id = bed$geneID,
  stringsAsFactors = FALSE
)

bin_bed_200k <- bin_bed_200k[bin_bed_200k[, 1] %in% paste0("chr", 1:22), ]

seg_filtered <- CreateSegtableNoWGS(
  mtx = Input_filtered$mtx,
  barcodes = Input_filtered$barcodes,
  features = Input_filtered$features,
  Cov_obj = Cov_obj,
  Zest = Zest,
  gtf = gtf,
  size = size,
  bin_bed = bin_bed_200k,
  celltype0 = celltype0,
  bin_mtx = NULL,
  plot_seg = TRUE,
  hmm_states = c(0.5, 1.5, 2),
  max_qt = 0.95,
  nmean = 500,
  rm_extreme = 2,
  adj = -0.5,
  dir_path = output_dir
)

second_seg_table <- data.frame(
  chr = seg_filtered[, 1],
  start = as.numeric(seg_filtered[, 2]),
  end = as.numeric(seg_filtered[, 3]),
  states = 1,
  length = as.numeric(seg_filtered[, 3]) - as.numeric(seg_filtered[, 2]),
  mean = 0, var = 0, Var1 = 1:nrow(seg_filtered), Freq = 50000,
  chrr = paste0(seg_filtered[, 1], ":", seg_filtered[, 2]),
  stringsAsFactors = FALSE
)

message(">>> [Clonalscope_NoWGS] Step 5: Running Round 2 (Refined)...")

Cov_obj2 <- RunCovCluster(
  mtx = Input_filtered$mtx,
  barcodes = Input_filtered$barcodes,
  features = Input_filtered$features,
  bed = bed,
  celltype0 = celltype0,
  var_pt = 0.99,
  var_pt_ctrl = 0.99,
  include = "all",
  alpha_source = "all",
  ctrl_region = NULL,
  seg_table_filtered = second_seg_table, 
  size = size,
  dir_path = output_dir,
  breaks = 50,
  ngene_filter=10,
  prep_mode = "intersect",
  est_cap = 2,
  clust_mode = "all",
  mincell = mincell
)



final_tn_res <- MalignantAssignment(Cov_obj2,cutoff=0.5)

final_assignment_full <- setNames(
  rep("Unknown", length(colnames(Input_filtered$mtx))),
  colnames(Input_filtered$mtx)
)

final_assignment_full[names(final_tn_res$final_assignment)] <-
  final_tn_res$final_assignment

final_assign_df <- data.frame(
  Barcode = names(final_assignment_full),
  Assignment = as.character(final_assignment_full),
  stringsAsFactors = FALSE
)

write.table(
  final_assign_df,
  file = file.path(output_dir, "clonalscope_final_tn_assignment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Save output
write.table(Cov_obj2$result_final$df_obj$df, 
            file = file.path(output_dir, "clonalscope_cell_seg_matrix.tsv"), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

zest_vec <- Cov_obj2$result_final$result$Zest
zest_df <- data.frame(Barcode = names(zest_vec), Cluster = as.numeric(zest_vec))
write.table(zest_df, 
            file = file.path(output_dir, "clonalscope_zest_clusters.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)   


annot_vec <- Cov_obj2$result_final$result$annot
# annot do not has barcode, use zest barcodes
annot_df <- data.frame(Barcode = names(zest_vec), Type = as.character(annot_vec))
write.table(annot_df, 
            file = file.path(output_dir, "clonalscope_cell_annotations.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)        

  
saveRDS(Cov_obj2, file.path(output_dir, "Cov_obj_second.rds"))
message(">>> [Clonalscope_NoWGS] Done.")