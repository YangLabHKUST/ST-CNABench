#!/usr/bin/env bash
set -euo pipefail

# Read positional arguments.
# Use ${1:-} so missing args fail explicitly later without tripping set -u here.
CALICOST_DIR="${1:-}"
CONFIG_YAML="${2:-}"
PURITY_CONFIG="${3:-}"
CNA_CONFIG="${4:-}"
PREPROC_DIR="${5:-}"
CORES="${6:-20}"
USE_TUMOR_PURITY_RAW="${7:-true}"

# Path
ABS_CALICOST_DIR=$(realpath "$CALICOST_DIR")
OUTPUT_DIR=$(dirname "$PREPROC_DIR")
PARSED_INPUT_DIR="${OUTPUT_DIR}/parsed_inputs"
TUMORPROP_FILE="${OUTPUT_DIR}/loh_estimator_tumor_prop.tsv"
HAS_PURITY=true

case "${USE_TUMOR_PURITY_RAW,,}" in
  true|1|yes|y)
    USE_TUMOR_PURITY=true
    ;;
  false|0|no|n)
    USE_TUMOR_PURITY=false
    ;;
  *)
    echo "[run_calicost.sh] WARNING: Invalid use_tumor_purity='$USE_TUMOR_PURITY_RAW'. Falling back to true."
    USE_TUMOR_PURITY=true
    ;;
esac

echo "[run_calicost.sh] use_tumor_purity = $USE_TUMOR_PURITY"

# --- Step 1: Preprocessing ---
if [ -f "$PREPROC_DIR/unique_snp_ids.npy" ]; then
  echo "[run_calicost.sh] Step 1 already done, skipping snakemake preprocessing."
else
  echo "[run_calicost.sh] Step 1: snakemake preprocessing..."
  # Workaround for snakemake source-cache path bug on some envs:
  # use a dataset-local tmp dir and a relative snakefile path.
  SNAKE_TMP_ROOT="${OUTPUT_DIR}/snakemake_tmp"
  mkdir -p "$SNAKE_TMP_ROOT"
  export TMPDIR="$SNAKE_TMP_ROOT"
  export SNAKEMAKE_TMPDIR="$SNAKE_TMP_ROOT"
  (
    cd "$ABS_CALICOST_DIR"
    snakemake \
      --cores "$CORES" \
      --configfile "$CONFIG_YAML" \
      --snakefile "calicost.smk" \
      --printshellcmds \
      --keep-going \
      all
  ) 2>&1
fi

# --- Step 2: Purity Estimation ---
if [ "$USE_TUMOR_PURITY" = "false" ]; then
  echo "[run_calicost.sh] Step 2: skipped (tumor purity disabled by config)."
  HAS_PURITY=false
  if [ -f "$TUMORPROP_FILE" ]; then
    echo "[run_calicost.sh] Removing stale $TUMORPROP_FILE because use_tumor_purity=false"
    rm -f "$TUMORPROP_FILE"
  fi
else
  if [ -f "$TUMORPROP_FILE" ]; then
    echo "[run_calicost.sh] Step 2 already done, skipping tumor purity estimation."
  else
    echo "[run_calicost.sh] Step 2: estimating tumor purity..."
    set +e
    OMP_NUM_THREADS=1 python "$ABS_CALICOST_DIR/src/calicost/estimate_tumor_proportion.py" -c "$PURITY_CONFIG"
    # FallBack to None if purity estimation fails. See: https://github.com/raphael-group/CalicoST/issues/18
    if [ $? -ne 0 ]; then
      echo "[run_calicost.sh] WARNING: Purity estimation failed. Falling back to None."
      HAS_PURITY=false
      rm -f "$TUMORPROP_FILE"
    fi
    set -e
  fi

  if [ ! -f "$TUMORPROP_FILE" ]; then
    HAS_PURITY=false
  fi
fi

if [ -d "$PARSED_INPUT_DIR" ]; then
  echo "[run_calicost.sh] Removing cached $PARSED_INPUT_DIR to force re-parse"
  rm -rf "$PARSED_INPUT_DIR"
else
  echo "[run_calicost.sh] No parsed_inputs found at $PARSED_INPUT_DIR, skipping removal"
fi

echo "[run_calicost.sh] Step 3: inferring CNAs..."
if [ "$HAS_PURITY" = "false" ]; then
  echo "[run_calicost.sh] Running Step 3 WITHOUT tumor purity file..."
  sed -i 's|^tumorprop_file.*|tumorprop_file : None|' "$CNA_CONFIG"
fi

OMP_NUM_THREADS=1 python "$ABS_CALICOST_DIR/src/calicost/calicost_main.py" -c "$CNA_CONFIG"
echo "[run_calicost.sh] Done."
