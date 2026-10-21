suppressPackageStartupMessages({
  library(dplyr); library(readr); library(pROC)
})

# Phase 8 build order, same as CHAR_ORDER in eval_build_train_input.py
CHARACTER_COLUMNS_P8 <- c(
  "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired", "AGly_acquired",
  "Col_acquired", "Tet_acquired", "MLS_acquired", "Rif_acquired", "Tmt_acquired",
  "Fcyn_acquired", "Phe_acquired", "Yersiniabactin", "Aerobactin", "Salmochelin",
  "Colibactin", "rmpA", "rmpA2", "Bla_acquired", "Tgc_acquired"
)
RES_CHARS <- c(
  "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired",
  "AGly_acquired", "Col_acquired", "Tet_acquired", "MLS_acquired",
  "Rif_acquired", "Tmt_acquired", "Fcyn_acquired", "Phe_acquired",
  "Bla_acquired", "Tgc_acquired"
)
VIR_CHARS <- c("Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2")

missing_set_for_state <- function(starting_state, char_order) {
  res_idx <- match(RES_CHARS, char_order) - 1L
  vir_idx <- match(VIR_CHARS, char_order) - 1L
  if (starting_state == "resistance_only") return(vir_idx)
  if (starting_state == "virulence_only")  return(res_idx)
  stop(sprintf("Unrecognized starting_state value: '%s'", starting_state))
}

read_routes <- function(path) {
  stopifnot(file.exists(path))
  as.list(readLines(path))
}

score_lineage <- function(current_state, missing_idx_set, routes) {
  present_idx <- which(current_state == 1L) - 1L
  hits <- 0L; valid <- 0L
  for (route in routes) {
    steps <- as.integer(strsplit(trimws(route), "\\s+")[[1]])
    pos <- match(present_idx, steps)
    if (any(is.na(pos))) next
    if (length(pos) == 0) next
    trunc_at <- max(pos)
    if (trunc_at >= length(steps)) next
    next_feat <- steps[trunc_at + 1]
    valid <- valid + 1L
    if (next_feat %in% missing_idx_set) hits <- hits + 1L
  }
  if (valid == 0) return(NA_real_)
  hits / valid
}

score_task_table <- function(task_tbl, routes, cm, char_order, train_cutoff = 2017) {
  task_tbl %>%
    rowwise() %>%
    mutate(
      predicted_risk = {
        row_state <- cm %>%
          filter(base_ST == .env$base_ST, collection_year <= .env$train_cutoff) %>%
          select(all_of(char_order)) %>%
          summarise(across(everything(), ~ as.integer(any(.x == 1L, na.rm = TRUE)))) %>%
          as.integer()
        missing_idx_set <- missing_set_for_state(starting_state, char_order)
        score_lineage(row_state, missing_idx_set, routes)
      }
    ) %>%
    ungroup()
}

report_auroc <- function(d, label) {
  scored <- d %>% filter(!is.na(predicted_risk))
  n_pos <- sum(scored$became_dual_risk == 1)
  n_neg <- sum(scored$became_dual_risk == 0)
  if (n_pos == 0 || n_neg == 0) {
    cat(sprintf("%s: cannot compute (pos=%d neg=%d)\n", label, n_pos, n_neg))
    return(invisible(NULL))
  }
  roc_obj <- roc(scored$became_dual_risk, scored$predicted_risk, quiet = TRUE, direction = "<")
  ci_val <- ci.auc(roc_obj, method = "delong")
  cat(sprintf("%s: AUROC=%.3f  95%% CI=[%.3f, %.3f]  n=%d (pos=%d)\n",
              label, as.numeric(auc(roc_obj)), ci_val[1], ci_val[3], nrow(scored), n_pos))
}

cm           <- read_csv("phase3_out/character_matrix.csv", show_col_types = FALSE)
validate_tbl <- read_csv("phase8_out/phase8_validate_task_table.csv", show_col_types = FALSE)
test_tbl     <- read_csv("phase8_out/phase8_test_task_table.csv", show_col_types = FALSE)

RUN_DIR <- "phase9_subsample_out"

for (i in 1:3) {
  route_path <- file.path(RUN_DIR, paste0("kp_subsample_rep", i, "-routes.txt"))
  if (!file.exists(route_path)) {
    cat(sprintf("Replicate %d: routes file not found yet (%s) — skipping\n", i, route_path))
    next
  }
  routes <- read_routes(route_path)
  cat(sprintf("\nReplicate %d: %d trajectories\n", i, length(routes)))

  val_scored  <- score_task_table(validate_tbl, routes, cm, CHARACTER_COLUMNS_P8)
  test_scored <- score_task_table(test_tbl, routes, cm, CHARACTER_COLUMNS_P8)

  write_csv(val_scored,  file.path(RUN_DIR, paste0("subsample_rep", i, "_validate_scored.csv")))
  write_csv(test_scored, file.path(RUN_DIR, paste0("subsample_rep", i, "_test_scored.csv")))

  report_auroc(val_scored,  sprintf("  Replicate %d VALIDATE", i))
  report_auroc(test_scored, sprintf("  Replicate %d TEST", i))
}

cat("\nBaseline for comparison — sensitivity-check (full TRAIN, 3-seed pooled): VALIDATE 0.669, TEST 0.685\n")