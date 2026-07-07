#!/bin/bash

# Arguments
SAMPLE_NAME=$1
BAM_FILE=$2
BARCODES_FILE=$3
OUTPUT_DIR=$4
N_CORES=$5
MATRIX_DIR=$6
META_FILE=$7
GENOME_VER=$8
PILEUP_SCRIPT=$9
EAGLE_PATH=${10}
GMAP_PATH=${11}
SNP_VCF=${12}
PANEL_DIR=${13}
UMI_TAG=${14}
CELL_TAG=${15}
N_CLONES=${16}
# cellsnp-lite output 
ALLELE_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_allele_counts.tsv.gz"

if [ -f "$ALLELE_FILE" ]; then
    echo ">>> [Numbat Wrapper] Allele file exists: $ALLELE_FILE. Skipping."
else
    echo ">>> [Numbat Wrapper] Running pileup_and_phase.R..."
    

    mkdir -p "$OUTPUT_DIR"

    Rscript "$PILEUP_SCRIPT" \
        --label "$SAMPLE_NAME" \
        --samples "$SAMPLE_NAME" \
        --bams "$BAM_FILE" \
        --barcodes "$BARCODES_FILE" \
        --outdir "$OUTPUT_DIR" \
        --gmap "$GMAP_PATH" \
        --snpvcf "$SNP_VCF" \
        --paneldir "$PANEL_DIR" \
        --ncores "$N_CORES" \
        --eagle "$EAGLE_PATH" \
        --UMItag "$UMI_TAG" \
        --cellTAG "$CELL_TAG"
        
    if [ ! -f "$ALLELE_FILE" ]; then
        echo "Error: Allele file not found at $ALLELE_FILE"
        exit 1
    fi
fi

echo ">>> [Numbat Wrapper] Running Numbat Analysis..."

SCRIPT_DIR=$(dirname "$0")

Rscript "$SCRIPT_DIR/_numbat_analysis.R" \
    "$MATRIX_DIR" \
    "$OUTPUT_DIR" \
    "$META_FILE" \
    "$ALLELE_FILE" \
    "$N_CLONES" \
    "$N_CORES" \
    "$GENOME_VER"

echo ">>> [Numbat Wrapper] Pipeline Finished."