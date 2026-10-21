#!/usr/bin/env Rscript
# Phase 9: geographic holdout split
#
# Decision rules fixed BEFORE looking at any downstream AUROC, using only the
# diagnostic output (country distribution, purity threshold):
#
#   HOLDOUT_COUNTRY = "China"
#     Single largest country by TRAIN-window genome count (21.5%, vs USA's 18% and
#     a long tail below that). An outcome-blind rule, not a pick among similar options.
#
#   PURITY_THRESHOLD = 0.8
#     Same cutoff the diagnostic script used to call a lineage "meaningfully mixed"
#     (240/729). Reused here, not chosen after seeing this script's results.
#
# The unit held out is the lineage, not the genome, to match
# robustness_lineage_holdout_split.R so the two robustness checks are comparable.
# A lineage is HELD_OUT only if its primary country is China AND that assignment is
# unambiguous (frac >= 0.8). Mixed lineages stay in TRAIN untouched.

set.seed(42)  # nothing random here (assignment is deterministic), kept for provenance

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

CHAR_MATRIX_PATH    <- "phase3_out/character_matrix.csv"
VALIDATE_TABLE_PATH <- "phase8_out/phase8_validate_task_table.csv"
TEST_TABLE_PATH     <- "phase8_out/phase8_test_task_table.csv"
OUT_DIR             <- "phase9_out"
dir.create(OUT_DIR, showWarnings = FALSE)

HOLDOUT_COUNTRY   <- "China"
PURITY_THRESHOLD  <- 0.8

cm           <- read_csv(CHAR_MATRIX_PATH, show_col_types = FALSE)
validate_tbl <- read_csv(VALIDATE_TABLE_PATH, show_col_types = FALSE)
test_tbl     <- read_csv(TEST_TABLE_PATH, show_col_types = FALSE)

stopifnot(all(c("base_ST", "country", "collection_year") %in% names(cm)))

train_genomes <- cm %>% filter(collection_year <= 2017)

# primary country per lineage; count ties are broken alphabetically by country
# so the result doesn't depend on row order
lineage_country <- train_genomes %>%
  count(base_ST, country) %>%
  group_by(base_ST) %>%
  mutate(total = sum(n), frac = n / total) %>%
  arrange(base_ST, desc(n), country) %>%
  slice_head(n = 1) %>%
  ungroup()

holdout_lineages <- lineage_country %>%
  filter(country == HOLDOUT_COUNTRY, frac >= PURITY_THRESHOLD) %>%
  pull(base_ST)

all_lineages   <- sort(unique(cm$base_ST))
train_lineages <- setdiff(all_lineages, holdout_lineages)

n_val_in_holdout  <- sum(validate_tbl$base_ST %in% holdout_lineages)
n_test_in_holdout <- sum(test_tbl$base_ST %in% holdout_lineages)

cat(sprintf("Holdout country: %s (purity threshold >= %.1f)\n", HOLDOUT_COUNTRY, PURITY_THRESHOLD))
cat(sprintf("Holdout lineages: %d / %d total\n", length(holdout_lineages), length(all_lineages)))
cat(sprintf("VALIDATE task-table lineages landing in HELD-OUT: %d / %d\n", n_val_in_holdout, nrow(validate_tbl)))
cat(sprintf("TEST task-table lineages landing in HELD-OUT: %d / %d\n", n_test_in_holdout, nrow(test_tbl)))

write_csv(tibble(base_ST = holdout_lineages, split = "HELD_OUT"),
          file.path(OUT_DIR, "geo_holdout_assignment.csv"))
write_csv(tibble(base_ST = train_lineages, split = "TRAIN_LINEAGE"),
          file.path(OUT_DIR, "geo_train_assignment.csv"))

writeLines(c(
  "PHASE 9 — Geographic holdout split, dated decision record",
  paste("Date:", Sys.Date()),
  paste("Holdout country:", HOLDOUT_COUNTRY),
  paste("Purity threshold:", PURITY_THRESHOLD),
  paste("Holdout lineages:", length(holdout_lineages), "/", length(all_lineages)),
  paste("VALIDATE task lineages in holdout:", n_val_in_holdout, "/", nrow(validate_tbl)),
  paste("TEST task lineages in holdout:", n_test_in_holdout, "/", nrow(test_tbl)),
  "",
  "Rationale: China selected as the single largest TRAIN-window country by",
  "genome count (21.5%), an objective outcome-blind criterion given no country",
  "or natural regional grouping captures a majority of the data. Purity",
  "threshold (0.8) reused from the diagnostic script's own reporting cutoff,",
  "not newly chosen after seeing this script's results. Mixed lineages",
  "(frac < 0.8) are left in TRAIN untouched, neither held out nor excluded."
), file.path(OUT_DIR, "phase9_geo_holdout_decision.txt"))

cat("\nWrote:\n")
cat(" -", file.path(OUT_DIR, "geo_holdout_assignment.csv"), "\n")
cat(" -", file.path(OUT_DIR, "geo_train_assignment.csv"), "\n")
cat(" -", file.path(OUT_DIR, "phase9_geo_holdout_decision.txt"), "\n")