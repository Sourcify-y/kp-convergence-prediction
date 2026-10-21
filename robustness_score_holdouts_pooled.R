#!/usr/bin/env Rscript
# Phase 9: 3-seed pooled route scoring
#
# Adapted from eval_score_lineages.R and adopted as the standard scoring method
# (decided before this run's results existed).
#
# Usage:  Rscript robustness_score_holdouts_pooled.R [holdout|geo|sensitivity-check]
#
#   holdout            scores the Phase 9 lineage-holdout run (s0/s1/s2 from
#                      phase9_out/kp_holdout_run_model1_s*) against the part of the Phase 8 task
#                      tables whose lineages are HELD-OUT. This is the primary Phase 9 result.
#   geo                same, for the geographic-holdout run (kp_geo_holdout_run_model1_s*)
#   sensitivity-check  re-scores the original Phase 8 model (phase8_out/kp_train_only_model1_s*)
#                      with 3-seed pooling against the full VALIDATE/TEST tables. Written to
#                      separate files (phase8_sensitivity_check_pooled_*.csv) and reported
#                      ALONGSIDE the original s0-only AUROCs, never in place of them.
#
# Pooling: all sampled trajectories from the 3 seeds' -routes.txt files are concatenated before
# truncating/tallying, rather than scoring each seed separately and averaging. Each seed's full
# route sample then counts as extra draws from the same posterior.
#
# Column order matters (fixed 2026-09): Phase 8's routes were built with
# eval_build_train_input.py's CHAR_ORDER, which is a different permutation from CHARACTER_COLUMNS
# below (the character_matrix.csv header order used for the Phase 9 holdout/geo runs). Scoring one
# with the other's order scrambled every index and gave ~80% spurious NAs. CHARACTER_COLUMNS_P8 is
# Phase 8's real order, and char_order is passed through score_task_table() and
# missing_set_for_state() so each mode uses the order that matches its routes' training data.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pROC)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else stop("Usage: Rscript robustness_score_holdouts_pooled.R [holdout|geo|sensitivity-check]")
if (!mode %in% c("holdout", "geo", "sensitivity-check")) {
  stop("mode must be 'holdout', 'geo', or 'sensitivity-check'")
}

# character_matrix.csv header order (confirmed 2026-09). Used for the Phase 9 holdout/geo runs,
# written and read consistently by robustness_*_holdout_train_input.py and this script.
CHARACTER_COLUMNS <- c(
  "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
  "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired",
  "Tet_acquired", "Tgc_acquired", "Tmt_acquired", "Bla_acquired",
  "Bla_ESBL_acquired", "Bla_Carb_acquired", "Yersiniabactin",
  "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"
)

RES_CHARS <- c(
  "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired",
  "AGly_acquired", "Col_acquired", "Tet_acquired", "MLS_acquired",
  "Rif_acquired", "Tmt_acquired", "Fcyn_acquired", "Phe_acquired",
  "Bla_acquired", "Tgc_acquired"
)
VIR_CHARS <- c(
  "Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"
)
stopifnot(length(RES_CHARS) == 14, length(VIR_CHARS) == 6,
          setequal(c(RES_CHARS, VIR_CHARS), CHARACTER_COLUMNS))

# Phase 8's real build order (confirmed 2026-09 against eval_build_train_input.py's CHAR_ORDER).
# Only for sensitivity-check mode. Do NOT use it for holdout/geo: those routes were built with
# CHARACTER_COLUMNS above.
CHARACTER_COLUMNS_P8 <- c(
  "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired", "AGly_acquired",
  "Col_acquired", "Tet_acquired", "MLS_acquired", "Rif_acquired", "Tmt_acquired",
  "Fcyn_acquired", "Phe_acquired", "Yersiniabactin", "Aerobactin", "Salmochelin",
  "Colibactin", "rmpA", "rmpA2", "Bla_acquired", "Tgc_acquired"
)
stopifnot(setequal(CHARACTER_COLUMNS_P8, CHARACTER_COLUMNS))  # same 20 characters, different order

# 0-based indices of the characters a lineage is still missing, given its starting state
missing_set_for_state <- function(starting_state, char_order) {
  res_idx <- match(RES_CHARS, char_order) - 1L
  vir_idx <- match(VIR_CHARS, char_order) - 1L
  if (starting_state == "resistance_only") return(vir_idx)
  if (starting_state == "virulence_only")  return(res_idx)
  stop(sprintf("Unrecognized starting_state value: '%s' (expected 'resistance_only' or 'virulence_only')", starting_state))
}

read_pooled_routes <- function(route_paths) {
  routes <- list()
  for (p in route_paths) {
    stopifnot(file.exists(p))
    seed_routes <- readLines(p)
    routes <- c(routes, as.list(seed_routes))
  }
  routes
}

score_lineage <- function(current_state, missing_idx_set, pooled_routes) {
  # Same method as eval_score_lineages.R: truncate each route at whichever of the lineage's
  # ALREADY-PRESENT characters comes latest, then look at what comes next. It doesn't require the
  # lineage's absent characters to be absent from the route up to that point. This is
  # intentionally the looser, original, pre-registered definition.
  present_idx <- which(current_state == 1L) - 1L  # 0-indexed
  hits <- 0L
  valid <- 0L
  for (route in pooled_routes) {
    steps <- as.integer(strsplit(trimws(route), "\\s+")[[1]])
    pos <- match(present_idx, steps)
    if (any(is.na(pos))) next
    trunc_at <- max(pos)
    if (trunc_at >= length(steps)) next
    next_feat <- steps[trunc_at + 1]
    valid <- valid + 1L
    if (next_feat %in% missing_idx_set) hits <- hits + 1L
  }
  if (valid == 0) return(NA_real_)
  hits / valid
}

score_task_table <- function(task_tbl, pooled_routes, cm, char_order, train_cutoff = 2017) {
  results <- task_tbl %>%
    rowwise() %>%
    mutate(
      predicted_risk = {
        row_state <- cm %>%
          filter(base_ST == .env$base_ST, collection_year <= .env$train_cutoff) %>%
          select(all_of(char_order)) %>%
          summarise(across(everything(), ~ as.integer(any(.x == 1L, na.rm = TRUE)))) %>%
          as.integer()
        missing_idx_set <- missing_set_for_state(starting_state, char_order)
        score_lineage(row_state, missing_idx_set, pooled_routes)
      }
    ) %>%
    ungroup()
  results
}

if (mode == "holdout") {
  cat("=== Mode: holdout (Phase 9 primary result) ===\n")
  HOLDOUT_RUN_DIR <- "phase9_out"
  holdout_assign <- read_csv("phase9_out/lineage_holdout_assignment.csv", show_col_types = FALSE)
  cm             <- read_csv("phase3_out/character_matrix.csv", show_col_types = FALSE)
  validate_tbl   <- read_csv("phase8_out/phase8_validate_task_table.csv", show_col_types = FALSE)
  test_tbl       <- read_csv("phase8_out/phase8_test_task_table.csv", show_col_types = FALSE)

  validate_holdout <- validate_tbl %>% filter(base_ST %in% holdout_assign$base_ST)
  test_holdout     <- test_tbl %>% filter(base_ST %in% holdout_assign$base_ST)

  cat(sprintf("Scoring %d held-out VALIDATE lineages, %d held-out TEST lineages\n",
              nrow(validate_holdout), nrow(test_holdout)))

  route_paths <- file.path(HOLDOUT_RUN_DIR, paste0("kp_holdout_run_model1_s", 0:2, "-routes.txt"))
  pooled <- read_pooled_routes(route_paths)
  cat(sprintf("Pooled %d trajectories across 3 seeds\n", length(pooled)))

  validate_scored <- score_task_table(validate_holdout, pooled, cm, CHARACTER_COLUMNS)
  test_scored     <- score_task_table(test_holdout, pooled, cm, CHARACTER_COLUMNS)

  dir.create("phase9_out", showWarnings = FALSE)
  write_csv(validate_scored, "phase9_out/phase9_validate_holdout_scored.csv")
  write_csv(test_scored, "phase9_out/phase9_test_holdout_scored.csv")

  auroc_val  <- auc(roc(validate_scored$became_dual_risk, validate_scored$predicted_risk, quiet = TRUE, direction = "<"))
  auroc_test <- auc(roc(test_scored$became_dual_risk, test_scored$predicted_risk, quiet = TRUE, direction = "<"))
  cat(sprintf("\nLineage-holdout VALIDATE AUROC: %.3f (n=%d)\n", auroc_val, nrow(validate_scored)))
  cat(sprintf("Lineage-holdout TEST AUROC:     %.3f (n=%d)\n", auroc_test, nrow(test_scored)))

} else if (mode == "geo") {
  cat("=== Mode: geo (Phase 9 geographic-holdout result) ===\n")
  GEO_RUN_DIR <- "phase9_out"
  geo_assign   <- read_csv("phase9_out/geo_holdout_assignment.csv", show_col_types = FALSE)
  cm           <- read_csv("phase3_out/character_matrix.csv", show_col_types = FALSE)
  validate_tbl <- read_csv("phase8_out/phase8_validate_task_table.csv", show_col_types = FALSE)
  test_tbl     <- read_csv("phase8_out/phase8_test_task_table.csv", show_col_types = FALSE)

  validate_geo <- validate_tbl %>% filter(base_ST %in% geo_assign$base_ST)
  test_geo     <- test_tbl %>% filter(base_ST %in% geo_assign$base_ST)

  cat(sprintf("Scoring %d held-out VALIDATE lineages, %d held-out TEST lineages\n",
              nrow(validate_geo), nrow(test_geo)))

  route_paths <- file.path(GEO_RUN_DIR, paste0("kp_geo_holdout_run_model1_s", 0:2, "-routes.txt"))
  pooled <- read_pooled_routes(route_paths)
  cat(sprintf("Pooled %d trajectories across 3 seeds\n", length(pooled)))

  validate_scored <- score_task_table(validate_geo, pooled, cm, CHARACTER_COLUMNS)
  test_scored     <- score_task_table(test_geo, pooled, cm, CHARACTER_COLUMNS)

  write_csv(validate_scored, "phase9_out/phase9_validate_geo_scored.csv")
  write_csv(test_scored, "phase9_out/phase9_test_geo_scored.csv")

  auroc_val  <- auc(roc(validate_scored$became_dual_risk, validate_scored$predicted_risk, quiet = TRUE, direction = "<"))
  auroc_test <- auc(roc(test_scored$became_dual_risk, test_scored$predicted_risk, quiet = TRUE, direction = "<"))
  cat(sprintf("\nGeo-holdout VALIDATE AUROC: %.3f (n=%d)\n", auroc_val, nrow(validate_scored)))
  cat(sprintf("Geo-holdout TEST AUROC:     %.3f (n=%d)\n", auroc_test, nrow(test_scored)))

} else {
  cat("=== Mode: sensitivity-check (Phase 8 re-scored, NOT a replacement) ===\n")
  cm           <- read_csv("phase3_out/character_matrix.csv", show_col_types = FALSE)
  validate_tbl <- read_csv("phase8_out/phase8_validate_task_table.csv", show_col_types = FALSE)
  test_tbl     <- read_csv("phase8_out/phase8_test_task_table.csv", show_col_types = FALSE)

  route_paths <- file.path("phase8_out", paste0("kp_train_only_model1_s", 0:2, "-routes.txt"))
  pooled <- read_pooled_routes(route_paths)
  cat(sprintf("Pooled %d trajectories across original Phase 8's 3 seeds\n", length(pooled)))

  validate_scored <- score_task_table(validate_tbl, pooled, cm, CHARACTER_COLUMNS_P8)
  test_scored     <- score_task_table(test_tbl, pooled, cm, CHARACTER_COLUMNS_P8)

  dir.create("phase8_out", showWarnings = FALSE)
  write_csv(validate_scored, "phase8_out/phase8_sensitivity_check_pooled_validate.csv")
  write_csv(test_scored, "phase8_out/phase8_sensitivity_check_pooled_test.csv")

  auroc_val  <- auc(roc(validate_scored$became_dual_risk, validate_scored$predicted_risk, quiet = TRUE, direction = "<"))
  auroc_test <- auc(roc(test_scored$became_dual_risk, test_scored$predicted_risk, quiet = TRUE, direction = "<"))
  cat(sprintf("\n3-seed-pooled VALIDATE AUROC: %.3f  (original s0-only: 0.664)\n", auroc_val))
  cat(sprintf("3-seed-pooled TEST AUROC:     %.3f  (original s0-only: 0.686)\n", auroc_test))
  cat("\nReport BOTH numbers in the manuscript, labeled: 's0-only (primary,\n")
  cat("pre-specified)' and '3-seed pooled (sensitivity check)'. Do not swap\n")
  cat("the primary number to whichever of these is larger.\n")
}