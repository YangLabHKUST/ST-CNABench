import argparse
import os
import sys
import pandas as pd
import shutil
from pathlib import Path
import logging

# environment variables to prevent err
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"

import xclone

def main():
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_name", required=True)
    parser.add_argument("--xcltk_dir", required=True) 
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--meta_file", required=True)
    parser.add_argument("--spatial_file", required=True)
    parser.add_argument("--n_clusters", required=True)
    args = parser.parse_args()

    # Paths
    xcltk_dir = Path(args.xcltk_dir)
    out_dir = Path(args.output_dir)
    ad_file = xcltk_dir / "3_baf_fc/xcltk.AD.mtx"
    dp_file = xcltk_dir / "3_baf_fc/xcltk.DP.mtx"
    baf_barcodes = xcltk_dir / "barcodes.tsv"

    n_clusters = int(args.n_clusters)
    rdr_mtx = xcltk_dir / "matrix.mtx"
    rdr_barcodes = xcltk_dir / "barcodes.tsv"
    rdr_features = xcltk_dir / "features.tsv"
    clean_anno_file = out_dir / "anno_file.clean.tsv"
    # Delete Comments and transfer "Normal" to "normal" "Tumor" to "tumor"
    with open(args.meta_file, 'r') as f_in, open(clean_anno_file, 'w') as f_out:
        for line in f_in:
            if not line.startswith('#'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) > 1:
                    for i in range(1, len(parts)):
                        parts[i] = parts[i].replace("Normal", "normal").replace("Tumor", "tumor")
                safe_line = '\t'.join(parts) + '\n'
                f_out.write(safe_line)

    logging.info(">>> [Xclone] Step 1: Loading Data...")
    # create data + extra annot + extra spatial
    # Ref: https://github.com/single-cell-genetics/XClone/issues/43
    BAF_adata = xclone.pp.xclonedata([str(ad_file), str(dp_file)], 'BAF',
                                     str(baf_barcodes),
                                     genome_mode="hg38_genes")
    
    BAF_adata = xclone.pp.extra_anno(BAF_adata, str(clean_anno_file), 
                                     barcodes_key="Barcode",
                                     cell_anno_key="tumor_normal", 
                                     sep="\t")
    
    BAF_adata = xclone.pp.extra_anno(BAF_adata, str(args.spatial_file), 
                                     barcodes_key="barcode", 
                                     cell_anno_key=None, sep=",")

    RDR_adata = xclone.pp.xclonedata(str(rdr_mtx), 'RDR',
                                     str(rdr_barcodes),
                                     genome_mode="hg38_genes")                             

    RDR_adata = xclone.pp.extra_anno(RDR_adata, str(clean_anno_file), 
                                     barcodes_key="Barcode",
                                     cell_anno_key="tumor_normal", 
                                     sep="\t")
    
    RDR_adata = xclone.pp.extra_anno(RDR_adata, str(args.spatial_file), 
                                     barcodes_key='barcode', 
                                     cell_anno_key=None, sep=",")


    run_out_dir = out_dir / "xclone_output"
    if not run_out_dir.exists(): run_out_dir.mkdir()
    
    logging.info(">>> [Xclone] Step 2: Running RDR Model...")
    # xclone config
    xconfig_rdr = xclone.XCloneConfig(dataset_name=args.sample_name, module="RDR", set_spatial=True)
    xconfig_rdr.outdir = str(run_out_dir)
    xconfig_rdr.cell_anno_key = "tumor_normal"
    xconfig_rdr.ref_celltype = "normal"
    xconfig_rdr.marker_group_anno_key = "tumor_normal"
    xconfig_rdr.top_n_marker = 15 ##
    xconfig_rdr.filter_ref_ave = 1.8 ##
    xconfig_rdr.min_gene_keep_num = 1000 ##
    xconfig_rdr.xclone_plot = True
    xconfig_rdr.plot_cell_anno_key = "tumor_normal"
    xconfig_rdr.fit_GLM_libratio = False
    xconfig_rdr.exclude_XY = False
    xconfig_rdr.remove_guide_XY = True
    
    RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file=xconfig_rdr)

    logging.info(">>> [Xclone] Step 3: Running BAF Model...")
    xconfig_baf = xclone.XCloneConfig(dataset_name=args.sample_name, module="BAF", set_spatial=True)
    xconfig_baf.update_info_from_rdr = False ## update using only BAF information.
    xconfig_baf.outdir = str(run_out_dir)
    xconfig_baf.cell_anno_key = "tumor_normal"
    xconfig_baf.ref_celltype = "normal"
    xconfig_baf.xclone_plot = True
    xconfig_baf.plot_cell_anno_key = "tumor_normal"
    xconfig_baf.phasing_region_key = "chr"
    xconfig_baf.update_info_from_rdr = False
    xconfig_baf.BAF_denoise = True
    
    BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, config_file=xconfig_baf)

    logging.info(">>> [Xclone] Step 4: Running Combine Model...")
    xconfig_cmb = xclone.XCloneConfig(dataset_name=args.sample_name, module="Combine", set_spatial=True)
    xconfig_cmb.outdir = str(run_out_dir)
    xconfig_cmb.cell_anno_key = "tumor_normal"
    xconfig_cmb.ref_celltype = "normal"
    xconfig_cmb.copygain_correct = False
    xconfig_cmb.xclone_plot = True
    xconfig_cmb.plot_cell_anno_key = "tumor_normal"
    xconfig_cmb.merge_loss = False
    xconfig_cmb.merge_loh = False
    xconfig_cmb.BAF_denoise = True
    xconfig_cmb.exclude_XY = True

    combine_Xdata = xclone.model.run_combine(RDR_Xdata,
                                             BAF_merge_Xdata,
                                             verbose=True,
                                             run_verbose=True,
                                             config_file=xconfig_cmb)
    
    # Detect normal and tumor based on BAF:
    logging.info(">>> [Xclone] Step 5: Analyzing tumor and normal...")
    xclone.al.exploreClustering(BAF_merge_Xdata, ref_anno_key = "tumor_normal", Xlayer = "posterior_mtx", max_clusters = 5)
    BAF_merge_Xdata_anno = xclone.al.OnestopBAFClustering(BAF_merge_Xdata, Xlayer = "posterior_mtx", 
                         n_clusters = 2, ref_anno_key = "tumor_normal", clone_anno_key = "clone(2)", 
                         plot_title = "XClone",
                         file_save_path = str(run_out_dir), file_save_name =  f"{args.sample_name}_2clones", complex_plot = False)

    # Subclone Detection based on combine data
    Comb_subclone_anno = xclone.al.OnestopBAFClustering(combine_Xdata, Xlayer = "prob1_merge", 
                         n_clusters = n_clusters, ref_anno_key = "tumor_normal", clone_anno_key = "subcluster", 
                         plot_title = "XClone",
                         file_save_path = str(run_out_dir), file_save_name =  f"{args.sample_name}_subcluster", complex_plot = False)
    logging.info(">>> [Xclone] Analysis Finished.")

if __name__ == "__main__":
    main()