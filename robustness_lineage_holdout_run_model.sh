#!/usr/bin/env bash
# Phase 9: lineage-holdout retrain.
# Run from the project root (or wherever hypertraps-ct/ is present).
#
# Flags, confirmed from the Phase 7 run banner:
#   --obs      observations file
#   --seed     random number seed
#   --length   length index
#   --kernel   kernel index (Phase 7/8 used 5)
#   --model    model (Phase 7/8 used 1, after --model 2 was ruled out as under-identified)
#   --label    output file prefix
# walkers/losses/APM/penalty are left at their defaults (200/0/0/0); Phase 7/8 didn't
# override them either.
#
# There is no --labels flag on the binary. kp-holdout-train-labels-20char.txt (written by
# robustness_lineage_holdout_train_input.py) is only a side file for mapping output rows
# back to genome ids later; it is not passed to hypertraps.ce.
set -euo pipefail

cd hypertraps-ct

OBS=Data/kp-holdout-train-data-20char.txt
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
    --label "${OUTDIR}/kp_holdout_run_model1_s${SEED}" \
    > "${OUTDIR}/run_log_holdout_s${SEED}.txt" 2>&1 &
  PIDS+=($!)
done
FAIL=0
for pid in "${PIDS[@]}"; do wait "$pid" || FAIL=1; done
[ "$FAIL" -eq 0 ] || { echo "ERROR: at least one seed failed; check the logs" >&2; exit 1; }

echo "All 3 seeds complete. Check seed agreement (model_compare_seeds.py) before trusting pooled routes."
echo "Logs: ${OUTDIR}/run_log_holdout_s{0,1,2}.txt"