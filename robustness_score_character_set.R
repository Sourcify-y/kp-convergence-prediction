suppressPackageStartupMessages({
  library(dplyr); library(readr); library(pROC)
})

RES_CHARS <- c(
  "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired",
  "AGly_acquired", "Col_acquired", "Tet_acquired", "MLS_acquired",
  "Rif_acquired", "Tmt_acquired", "Fcyn_acquired", "Phe_acquired",
  "Bla_acquired", "Tgc_acquired"
)
VIR_CHARS <- c("Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2")
COARSE_CHAR_ORDER <- c(RES_CHARS, "Hypervirulence_any")  # 15 columns, same as the build script

# 0-based indices of the characters a lineage is still missing, given its starting state
target_idx_for_state <- function(starting_state, char_order) {
  res_idx <- match(RES_CHARS, char_order) - 1L
  hv_idx  <- match("Hypervirulence_any", char_order) - 1L
  if (starting_state == "resistance_only") return(hv_idx)
  if (starting_state == "virulence_only")  return(res_idx)
  stop(sprintf("Unrecognized starting_state value: '%s'", starting_state))
}

read_routes <- function(path) { stopifnot(file.exists(path)); as.list(readLines(path)) }

# fraction of usable routes where the next trait after everything the lineage
# already has is in the target set
score_lineage <- function(current_state, target_idx, routes) {
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
    if (next_feat %in% target_idx) hits <- hits + 1L
  }
  if (valid == 0) return(NA_real_)
  hits / valid
}

score_task_table_coarse <- function(task_tbl, routes, cm, train_cutoff = 2017) {
  task_tbl %>%
    rowwise() %>%
    mutate(
      predicted_risk = {
        row_state <- cm %>%
          filter(base_ST == .env$base_ST, collection_year <= .env$train_cutoff) %>%
          select(all_of(COARSE_CHAR_ORDER)) %>%
          summarise(across(everything(), ~ as.integer(any(.x == 1L, na.rm = TRUE)))) %>%
          as.integer()
        target_idx <- target_idx_for_state(starting_state, COARSE_CHAR_ORDER)
        score_lineage(row_state, target_idx, routes)
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

# collapse the six virulence loci into one "any hypervirulence marker" column
cm <- cm %>%
  mutate(Hypervirulence_any = as.integer(rowSums(across(all_of(VIR_CHARS)) == 1L, na.rm = TRUE) > 0))

RUN_DIR <- "phase9_coarse_out"
pooled_routes <- list()
for (s in 0:2) {
  route_path <- file.path(RUN_DIR, paste0("kp_coarse_s", s, "-routes.txt"))
  if (!file.exists(route_path)) {
    cat(sprintf("Seed %d routes not found yet (%s) — skipping\n", s, route_path))
    next
  }
  pooled_routes <- c(pooled_routes, read_routes(route_path))
}
cat(sprintf("Pooled %d trajectories\n", length(pooled_routes)))

if (length(pooled_routes) > 0) {
  val_scored  <- score_task_table_coarse(validate_tbl, pooled_routes, cm)
  test_scored <- score_task_table_coarse(test_tbl, pooled_routes, cm)

  write_csv(val_scored,  file.path(RUN_DIR, "coarse_validate_scored.csv"))
  write_csv(test_scored, file.path(RUN_DIR, "coarse_test_scored.csv"))

  report_auroc(val_scored,  "Coarse (15-char) VALIDATE")
  report_auroc(test_scored, "Coarse (15-char) TEST")
}

cat("\nBaseline for comparison — 20-char sensitivity-check (full TRAIN, 3-seed pooled): VALIDATE 0.669, TEST 0.685\n")