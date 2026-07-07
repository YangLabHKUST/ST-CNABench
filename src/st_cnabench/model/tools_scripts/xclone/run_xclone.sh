#!/bin/bash

# Arguments
SAMPLE_NAME=$1
BAM_FILE=$2
BARCODES_FILE=$3
OUTPUT_DIR=$4
N_CORES=$5
MATRIX_DIR=$6
META_FILE=$7
SPATIAL_FILE=$8
SNP_VCF=$9
GENE_REGION=${10}
GMAP_PATH=${11}
EAGLE_PATH=${12}
PANEL_DIR=${13}
UMI_TAG=${14}
CELL_TAG=${15}
minCOUNT=${16}
minMAF=${17}
N_CLUSTERS=${18}
# Output Dir
XCLTK_DIR="$OUTPUT_DIR/xcltk_output"
mkdir -p "$XCLTK_DIR"

BAF_OUT_MARKER="$XCLTK_DIR/3_baf_fc/xcltk.AD.mtx"

# Check if exists
if [ -f "$BAF_OUT_MARKER" ]; then
    echo ">>> [Xclone Wrapper] xcltk BAF output exists. Skipping..."
else
    echo ">>> [Xclone Wrapper] Running xcltk baf..."
    
    xcltk baf \
        --label "$SAMPLE_NAME" \
        --sam "$BAM_FILE" \
        --barcode "$BARCODES_FILE" \
        --snpvcf "$SNP_VCF" \
        --region "$GENE_REGION" \
        --outdir "$XCLTK_DIR" \
        --gmap "$GMAP_PATH" \
        --eagle "$EAGLE_PATH" \
        --paneldir "$PANEL_DIR" \
        --ncores "$N_CORES" \
        --cellTAG "$CELL_TAG" \
        --UMItag "$UMI_TAG" \
        --minCOUNT "$minCOUNT" \
        --minMAF "$minMAF" 

    if [ ! -f "$BAF_OUT_MARKER" ]; then
        echo "Error: xcltk baf failed. Output not found."
        exit 1
    fi
fi

RDR_OUT_MARKER="$XCLTK_DIR/matrix.mtx"

if [ -f "$RDR_OUT_MARKER" ]; then
    echo ">>> [Xclone Wrapper] xcltk basefc output exists. Skipping..."
else
    echo ">>> [Xclone Wrapper] Running xcltk basefc..."
    
    xcltk basefc \
        --sam "$BAM_FILE" \
        --barcode "$BARCODES_FILE" \
        --region "$GENE_REGION" \
        --outdir "$XCLTK_DIR" \
        --ncores "$N_CORES" \
        --cellTAG "$CELL_TAG" \
        --UMItag "$UMI_TAG"
        
    if [ ! -f "$RDR_OUT_MARKER" ]; then
        echo "Error: xcltk basefc failed."
        exit 1
    fi
fi

echo ">>> [Xclone Wrapper] Preprocessing Done."

echo ">>> [Xclone Wrapper] Running Xclone Analysis..."

SCRIPT_DIR=$(dirname "$0")

python "$SCRIPT_DIR/_xclone_analysis.py" \
    --sample_name "$SAMPLE_NAME" \
    --xcltk_dir "$XCLTK_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --meta_file "$META_FILE" \
    --spatial_file "$SPATIAL_FILE" \
    --n_clusters "$N_CLUSTERS"

echo ">>> [Xclone Wrapper] Pipeline Finished."