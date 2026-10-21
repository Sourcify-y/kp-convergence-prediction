#!/usr/bin/env bash
# Phase 9: geographic-holdout retrain, same flags/config as the lineage holdout run.
# Runs seeds 0, 1, 2 in parallel.
set -euo pipefail

cd hypertraps-ct

OBS=Data/kp-geo-holdout-train-data-20char.txt
OUTDIR=../phase9_out
mkdir -p "$OUTDIR"

PIDS=()
for SEED in 0 1 2; do
  echo "=== Running seed $SEED ==="
  ./hypertraps.ce \
    --obs "$OBS" \
    --seed "$SEED" \
    --length 4 \
    --kernel 5 \
    --model 1 \
    --label "${OUTDIR}/kp_geo_holdout_run_model1_s${SEED}" \
    > "${OUTDIR}/run_log_geo_holdout_s${SEED}.txt" 2>&1 &
  PIDS+=($!)
done
FAIL=0
for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
[ "$FAIL" -eq 0 ] || { echo "ERROR: at least one seed failed; check the logs" >&2; exit 1; }

echo "All 3 seeds complete. Check seed agreement (model_compare_seeds.py) before trusting pooled routes."
echo "Logs: ${OUTDIR}/run_log_geo_holdout_s{0,1,2}.txt"