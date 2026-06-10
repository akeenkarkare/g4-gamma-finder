#!/bin/bash
# MC-init U-Net validation of the crystal-count knee (S=4 vs pP=5 vs h27=6).
# Confirms the knee ranking holds under the paper-grade model, not just the MLP.
#
# Run from the project root, on charger. `caffeinate` keeps the Mac awake so it
# can run with the lid closed (must be plugged in for clamshell mode).
#
#   ./run_unet_validation.sh
#
# Tune cost via env: SEEDS=5 EPOCHS=300 ./run_unet_validation.sh
set -e
cd "$(dirname "$0")"
source ml_env/bin/activate

run() {  # $1=label $2=dataset $3=filters
  echo "=================================================="
  caffeinate -s python3 train_unet.py --data "$2" --filters "$3" --label "$1"
}

echo "MC-init U-Net knee validation (SEEDS=${SEEDS:-3}, EPOCHS=${EPOCHS:-200})"
run "S (4)"   datasets_keep/dsN_S_nt_pixels.csv      build/filters_S.npy
run "pP (5)"  datasets_keep/dsShape_pP_nt_pixels.csv build/filters_pP.npy
run "h27 (6)" datasets_keep/dsHex_h27_nt_pixels.csv  build/filters_h27.npy
echo "=================================================="
echo "Done. Knee holds if  pP < h27 < S  (lower deg = better)."
