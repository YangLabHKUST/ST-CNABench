#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
  stop("Usage: Rscript run_clonalscope_wgs.R <matrix_dir> <output_dir> <meta_file> <wgs_tumor> <wgs_normal> <aux_dir> <gene_bed> <sample_name> <hmm_states> <mincell>")
}

matrix_dir   <- args[1]
output_dir   <- args[2]
meta_file    <- args[3]
wgs_tumor_f  <- args[4]
wgs_normal_f <- args[5]
aux_dir      <- args[6]  
gene_bed     <- args[7]
sample_name  <- args[8]
hmm_str      <- args[9]
mincell       <- args[10]

hmm_states <- as.numeric(unlist(strsplit(hmm_str, ",")))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

suppressPackageStartupMessages({
  library(Clonalscope)
  library(Matrix)
  library(data.table)
})

message(">>> [Clonalscope_WGS] Step 1: Loading Data...")

size <- read.table(file.path(aux_dir, "sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = FALSE)
cyclegenes <- readRDS(file.path(aux_dir, "cyclegenes.rds"))

bed <- read.table(gene_bed, sep='\t', header = TRUE)

if(ncol(bed) == 4) {
    bed <- bed[, c(2, 3, 4, 1)]
    colnames(bed) <- c("chr", "start", "end", "geneID")
} else {
    stop("Gene annotation file format error.")
}

#read bedgraph helper function
read_bedg_for_clonalscope <- function(path, add_start_plus1 = TRUE) {
  # bedgraph: chr, start, end, count (no header)
  dt <- fread(path, header = FALSE)
  setnames(dt, c("chr", "start", "end", "count"))

  dt[, chr := as.character(chr)]
  dt[, chr := sub("^chr", "", chr, ignore.case = TRUE)]

  # exclude X, Y, M
  dt <- dt[!(chr %in% c("X", "Y", "M", "MT"))]
  dt <- dt[chr %in% as.character(1:22)]

  dt[, chr_int := as.integer(chr)]
  setorder(dt, chr_int, start)

  # Start+1
  if (add_start_plus1) {
    dt[, start1 := as.integer(start) + 1L]
  } else {
    dt[, start1 := as.integer(start)]
  }
  dt[, bin_id := paste0("chr", chr, "-", start1, "-", as.integer(end))]
  out <- data.frame(count = dt$count)
  rownames(out) <- dt$bin_id
  return(out)
}

message(">>> [Clonalscope_WGS] Step 2: Processing WGS data...")
if (!file.exists(wgs_tumor_f) || !file.exists(wgs_normal_f)) {
    stop("WGS tumor or normal bedgraph files not found.")
}

WGSt <- read_bedg_for_clonalscope(wgs_tumor_f,  add_start_plus1 = TRUE)
WGSn <- read_bedg_for_clonalscope(wgs_normal_f, add_start_plus1 = TRUE)

# Intersect bins
common_bins <- intersect(rownames(WGSt), rownames(WGSn))
WGSt <- WGSt[common_bins, , drop = FALSE]
WGSn <- WGSn[common_bins, , drop = FALSE]
if(nrow(WGSt) == 0) stop("No common bins between Tumor and Normal WGS data!")
stopifnot(identical(rownames(WGSt), rownames(WGSn)))

message(sprintf(">>> [Clonalscope_WGS] Loaded %d bins.", nrow(WGSt)))

Obj_filtered <- Createobj_bulk(
    raw_counts = WGSt,
    ref_counts = WGSn, 
    samplename = sample_name,
    genome_assembly = "GRCh38", 
    dir_path = output_dir, 
    size = size, 
    assay = 'WGS'
)

# Segmentation using WGS/WES
Obj_filtered <- Segmentation_bulk(
    Obj_filtered = Obj_filtered,
    plot_seg = FALSE,
    hmm_states = hmm_states
)

message(">>> [Clonalscope_WGS] Step 3: Loading 10x Matrix...")
mtx_path <- file.path(matrix_dir, "matrix.mtx.gz")
bar_path <- file.path(matrix_dir, "barcodes.tsv.gz")
feat_path <- file.path(matrix_dir, "features.tsv.gz")

mtx <- readMM(mtx_path)
barcodes <- read.table(bar_path, stringsAsFactors = FALSE, sep='\t', header=FALSE)
features <- read.table(feat_path, stringsAsFactors = FALSE, sep='\t', header=FALSE)
gene_symbols = features[,2]
colnames(barcodes) <- "barcode"
colnames(features) <- c("gene_id","gene_name","feature_type")

rownames(mtx) <- features$gene_id   
colnames(mtx) <- barcodes$barcode
# Annot
if (is.null(meta_file) || meta_file == "None") {
  message(">>> [Clonalscope_WGS] No meta_file provided, will infer normal cells/spots from data.")
  initial_normal_spots = FindNormalReference(counts=mtx, gene_symbols=gene_symbols,method="pca")
  celltype0 = barcodes; rownames(celltype0) = barcodes[,1]
  celltype0$type = "tumor"
  celltype0[initial_normal_spots,"type"] = "normal"
  message(sprintf(">>> [Clonalscope_WGS] Inferred %d normal cells/spots.", length(initial_normal_spots)))
} else {
  message(">>> [Clonalscope_WGS] Meta_file provided, loading cell type annotation...")
  celltype0 <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(celltype0) <- c("barcode", "type") 
  celltype0$type <- ifelse(grepl("Tumor|tumor", celltype0$type, ignore.case = TRUE), "tumor", "normal")


  # Intersect cells
  common_cells <- intersect(colnames(mtx), celltype0$barcode)
  message(sprintf(">>> [Data] Common cells: %d", length(common_cells)))
  if (length(common_cells) == 0) stop("No common cells between matrix and metadata!")

  mtx <- mtx[, common_cells, drop = FALSE]
  barcodes <- data.frame(barcode = colnames(mtx), stringsAsFactors = FALSE)
  celltype0 <- celltype0[celltype0$barcode %in% common_cells, ]
}
#clustering_barcodes <- celltype0[celltype0$type == "tumor", 1]

#message(sprintf(">>> [Data] Tumor cells for clustering: %d", length(clustering_barcodes)))
#print(celltype0)

Input_filtered <- FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

message(">>> [Clonalscope_WGS] Step 4: Running RunCovCluster...")


Cov_obj <- RunCovCluster(
  mtx = Input_filtered$mtx, 
  barcodes = Input_filtered$barcodes,
  features = Input_filtered$features, 
  bed = bed,
  celltype0 = celltype0,
  var_pt = 0.99, 
  var_pt_ctrl = 0.99, 
  include = 'all',
  alpha_source = 'all', 
  ctrl_region = NULL,
  seg_table_filtered = Obj_filtered$seg_table,
  size = size, 
  dir_path = output_dir, 
  breaks = 50,
  ngene_filter=10,
  prep_mode = 'intersect', 
  est_cap = 2,
  clust_mode = 'cna_only',
  mincell = mincell
)

final_tn_res <- MalignantAssignment(Cov_obj,cutoff=0.5)

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

write.table(Cov_obj$result_final$df_obj$df, 
            file = file.path(output_dir, "clonalscope_cell_seg_matrix.tsv"), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

zest_vec <- Cov_obj$result_final$result$Zest
zest_df <- data.frame(Barcode = names(zest_vec), Cluster = as.numeric(zest_vec))
write.table(zest_df, 
            file = file.path(output_dir, "clonalscope_zest_clusters.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)    

annot_vec <- Cov_obj$result_final$result$annot
# annot do not has barcode, use zest barcodes
annot_df <- data.frame(Barcode = names(zest_vec), Type = as.character(annot_vec))
write.table(annot_df, 
            file = file.path(output_dir, "clonalscope_cell_annotations.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)        
saveRDS(Cov_obj, file.path(output_dir, "clonalscope_wgs_obj.rds"))
message(">>> [Clonalscope_WGS] Done.")