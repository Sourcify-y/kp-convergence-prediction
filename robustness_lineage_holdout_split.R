#!/usr/bin/env Rscript
# Phase 9: lineage holdout split
#
# Design (decided before this run existed, same endpoint-locking approach as the
# rest of the project): stratified 20% lineage holdout, weighted so the 143-lineage
# VALIDATE task table (and the 81-lineage TEST table) both keep enough held-out
# members to score against, without rebuilding the outcome tables from scratch.
#
# Expects a base_ST column in character_matrix.csv, phase8_validate_task_table.csv
# and phase8_test_task_table.csv.
#
# Target: ~40-45 of the 143 VALIDATE task-table lineages in HELD-OUT. Task-table
# lineages are oversampled relative to their share of the full 1,142-lineage pool
# (143/1142 ≈ 12.5%), because a flat 20% draw would only be expected to catch
# ~12.5% * 228 ≈ 28-29 of them, too few for a usable holdout task table. This is a
# stratification choice, not a result-driven one: the target count (40-45) was fixed
# before drawing, independent of any performance number.

set.seed(42)  # same seed as the Phase 8 paired bootstrap

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

CHAR_MATRIX_PATH    <- "phase3_out/character_matrix.csv"
VALIDATE_TABLE_PATH <- "phase8_out/phase8_validate_task_table.csv"
TEST_TABLE_PATH     <- "phase8_out/phase8_test_task_table.csv"
OUT_DIR             <- "phase9_out"
dir.create(OUT_DIR, showWarnings = FALSE)

cm           <- read_csv(CHAR_MATRIX_PATH, show_col_types = FALSE)
validate_tbl <- read_csv(VALIDATE_TABLE_PATH, show_col_types = FALSE)
test_tbl     <- read_csv(TEST_TABLE_PATH, show_col_types = FALSE)

stopifnot("base_ST" %in% names(cm))
stopifnot("base_ST" %in% names(validate_tbl))
stopifnot("base_ST" %in% names(test_tbl))

all_lineages     <- sort(unique(cm$base_ST))
n_lineages       <- length(all_lineages)
task_lineages    <- sort(unique(validate_tbl$base_ST))
n_task           <- length(task_lineages)
nontask_lineages <- setdiff(all_lineages, task_lineages)

cat(sprintf("Total lineages: %d\n", n_lineages))
cat(sprintf("Task-table (VALIDATE) lineages: %d\n", n_task))

# targets
TARGET_TOTAL_FRAC <- 0.20
TARGET_TASK_MIN   <- 40
TARGET_TASK_MAX   <- 45
target_task_n  <- sample(TARGET_TASK_MIN:TARGET_TASK_MAX, 1)  # fixed by range, not tuned post-hoc
target_total_n <- round(TARGET_TOTAL_FRAC * n_lineages)

if (target_task_n > n_task) {
  stop("target_task_n exceeds available task lineages — check TARGET_TASK_MIN/MAX against n_task")
}

# draw (keep this call order, it determines the seeded result)
holdout_task    <- sample(task_lineages, target_task_n)
remaining_slots <- max(target_total_n - target_task_n, 0)
holdout_nontask <- sample(nontask_lineages, min(remaining_slots, length(nontask_lineages)))

holdout_lineages <- sort(c(holdout_task, holdout_nontask))
train_lineages   <- setdiff(all_lineages, holdout_lineages)

# sanity checks (printed, not silently assumed)
n_val_in_holdout  <- sum(validate_tbl$base_ST %in% holdout_lineages)
n_test_in_holdout <- sum(test_tbl$base_ST %in% holdout_lineages)

cat(sprintf("Holdout lineages drawn: %d (%.1f%% of %d)\n",
            length(holdout_lineages), 100 * length(holdout_lineages) / n_lineages, n_lineages))
cat(sprintf("VALIDATE task-table lineages landing in HELD-OUT: %d / %d\n", n_val_in_holdout, nrow(validate_tbl)))
cat(sprintf("TEST task-table lineages landing in HELD-OUT: %d / %d\n", n_test_in_holdout, nrow(test_tbl)))

if (n_val_in_holdout < TARGET_TASK_MIN || n_val_in_holdout > TARGET_TASK_MAX) {
  warning("VALIDATE holdout count fell outside the 40-45 target band — re-check draw before proceeding.")
}

# outputs
write_csv(tibble(base_ST = holdout_lineages, split = "HELD_OUT"),
          file.path(OUT_DIR, "lineage_holdout_assignment.csv"))
write_csv(tibble(base_ST = train_lineages, split = "TRAIN_LINEAGE"),
          file.path(OUT_DIR, "lineage_train_assignment.csv"))

writeLines(c(
  "PHASE 9 — Lineage holdout split, dated decision record",
  paste("Date:", Sys.Date()),
  paste("Seed: 42"),
  paste("Total lineages:", n_lineages),
  paste("Target holdout fraction: 20%"),
  paste("Target task-table lineages in holdout: 40-45 (drawn:", target_task_n, ")"),
  paste("Actual holdout size:", length(holdout_lineages)),
  paste("VALIDATE task lineages in holdout:", n_val_in_holdout, "/", nrow(validate_tbl)),
  paste("TEST task lineages in holdout:", n_test_in_holdout, "/", nrow(test_tbl)),
  "",
  "Rationale: task-table lineages oversampled relative to base rate (143/1142)",
  "so the holdout retains a usable-sized scoring set without rebuilding outcome",
  "tables. Target band (40-45) and total fraction (20%) were fixed before the",
  "random draw; the draw itself (seed=42) was not repeated or cherry-picked."
), file.path(OUT_DIR, "phase9_lineage_holdout_decision.txt"))

cat("\nWrote:\n")
cat(" -", file.path(OUT_DIR, "lineage_holdout_assignment.csv"), "\n")
cat(" -", file.path(OUT_DIR, "lineage_train_assignment.csv"), "\n")
cat(" -", file.path(OUT_DIR, "phase9_lineage_holdout_decision.txt"), "\n")