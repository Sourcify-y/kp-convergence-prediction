#!/usr/bin/env Rscript
# Phase 9: geographic distribution diagnostic (run BEFORE designing the split)
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

cm           <- read_csv("phase3_out/character_matrix.csv", show_col_types = FALSE)
validate_tbl <- read_csv("phase8_out/phase8_validate_task_table.csv", show_col_types = FALSE)
test_tbl     <- read_csv("phase8_out/phase8_test_task_table.csv", show_col_types = FALSE)

stopifnot("country" %in% names(cm))

# genome-level country counts, TRAIN window only (that's what would be kept or excluded)
train_genomes <- cm %>% filter(collection_year <= 2017)
cat("=== Genome-level country distribution, TRAIN window (<=2017) ===\n")
train_genomes %>% count(country, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print(n = 30)

# each base_ST gets a "primary country" = plurality of its TRAIN genomes.
# Lineages can be split across countries, so also check how pure that assignment is
lineage_country <- train_genomes %>%
  count(base_ST, country) %>%
  group_by(base_ST) %>%
  mutate(total = sum(n), frac = n / total) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\n=== Lineage-level primary-country distribution ===\n")
lineage_country %>% count(country, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% print(n = 30)

cat("\n=== How 'pure' is the primary-country assignment? ===\n")
cat(sprintf("Lineages that are 100%% one country: %d / %d\n",
            sum(lineage_country$frac == 1), nrow(lineage_country)))
cat(sprintf("Lineages with primary-country fraction < 0.8 (meaningfully mixed): %d\n",
            sum(lineage_country$frac < 0.8)))

# overlay on the task-table lineages, since those are what has to be scoreable in the holdout region
cat("\n=== Task-table lineages by primary country ===\n")
val_country <- validate_tbl %>% left_join(lineage_country %>% select(base_ST, country), by = "base_ST")
test_country <- test_tbl %>% left_join(lineage_country %>% select(base_ST, country), by = "base_ST")
cat("VALIDATE task-table:\n")
val_country %>% count(country, sort = TRUE) %>% print(n = 20)
cat("TEST task-table:\n")
test_country %>% count(country, sort = TRUE) %>% print(n = 20)