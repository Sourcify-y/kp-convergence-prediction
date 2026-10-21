# eval_score_lineages.R
library(data.table)

CHAR_ORDER <- c("Bla_ESBL_acquired","Bla_Carb_acquired","Flq_acquired","Sul_acquired","AGly_acquired",
                "Col_acquired","Tet_acquired","MLS_acquired","Rif_acquired","Tmt_acquired",
                "Fcyn_acquired","Phe_acquired","Yersiniabactin","Aerobactin","Salmochelin",
                "Colibactin","rmpA","rmpA2","Bla_acquired","Tgc_acquired")

res_chars <- c("AGly_acquired","Bla_acquired","Bla_ESBL_acquired","Bla_Carb_acquired",
               "Col_acquired","Fcyn_acquired","Flq_acquired","MLS_acquired","Phe_acquired",
               "Rif_acquired","Sul_acquired","Tet_acquired","Tgc_acquired","Tmt_acquired")
vir_chars <- c("Yersiniabactin","Aerobactin","Salmochelin","Colibactin","rmpA","rmpA2")

# 0-indexed to match the routes file
res_idx <- which(CHAR_ORDER %in% res_chars) - 1
vir_idx <- which(CHAR_ORDER %in% vir_chars) - 1

# simulated routes from s0: 790 rows (79 posterior samples x 10), 20 cols, one full ordering per row
routes <- as.matrix(fread("phase8_out/kp_train_only_model1_s0-routes.txt", header = FALSE))
cat("Routes value range:", range(routes), "- dims:", dim(routes), "\n")
if (min(routes) == 1) {
  cat("Detected 1-indexed routes, converting to 0-indexed\n")
  routes <- routes - 1
}

# each lineage's TRAIN-time presence/absence over the 20 characters
dat <- fread("phase3_out/character_matrix.csv")
train_dat <- dat[collection_year <= 2017]
st_present <- train_dat[, lapply(.SD, function(x) as.integer(any(x == 1, na.rm = TRUE))),
                        by = base_ST, .SDcols = CHAR_ORDER]

# fraction of usable routes where the next feature acquired (after everything
# the lineage already has) is in the target set
score_st <- function(present_idx, target_idx) {
  hits <- 0; valid <- 0
  for (r in seq_len(nrow(routes))) {
    route <- routes[r, ]
    pos <- match(present_idx, route)
    if (any(is.na(pos))) next
    trunc_at <- max(pos)
    if (trunc_at >= length(route)) next
    next_feat <- route[trunc_at + 1]
    valid <- valid + 1
    if (next_feat %in% target_idx) hits <- hits + 1
  }
  if (valid == 0) return(NA_real_)
  hits / valid
}

score_task <- function(task_path, out_path) {
  task <- fread(task_path)
  task[, predicted_risk := NA_real_]
  for (i in seq_len(nrow(task))) {
    lin <- st_present[base_ST == task$base_ST[i]]
    if (nrow(lin) == 0) next
    present_idx <- which(as.integer(lin[1, ..CHAR_ORDER]) == 1) - 1
    target_idx <- if (task$starting_state[i] == "virulence_only") res_idx else vir_idx
    task$predicted_risk[i] <- score_st(present_idx, target_idx)
  }
  fwrite(task, out_path)
  cat(sprintf("%s: scored %d/%d lineages\n", out_path, sum(!is.na(task$predicted_risk)), nrow(task)))
}

score_task("phase8_out/phase8_validate_task_table.csv", "phase8_out/phase8_validate_scored.csv")
score_task("phase8_out/phase8_test_task_table.csv", "phase8_out/phase8_test_scored.csv")

# rank-based AUROC
calc_auc <- function(labels, scores) {
  ok <- !is.na(scores); labels <- labels[ok]; scores <- scores[ok]
  n1 <- sum(labels == 1); n0 <- sum(labels == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  r <- rank(scores)
  (sum(r[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

for (win in c("validate", "test")) {
  d <- fread(sprintf("phase8_out/phase8_%s_scored.csv", win))
  cat(sprintf("%s AUROC = %.3f (Phase 6 baseline: 0.590)\n", win, calc_auc(d$became_dual_risk, d$predicted_risk)))
}

calc_auprc <- function(labels, scores) {
  ok <- !is.na(scores); labels <- labels[ok]; scores <- scores[ok]
  o <- order(scores, decreasing = TRUE)
  labels <- labels[o]
  tp <- cumsum(labels); fp <- cumsum(1 - labels)
  precision <- tp / (tp + fp); recall <- tp / sum(labels)
  sum(diff(c(0, recall)) * precision)
}

boot_ci <- function(labels, scores, metric_fn, n_boot = 2000, seed = 42) {
  set.seed(seed)
  ok <- !is.na(scores); labels <- labels[ok]; scores <- scores[ok]
  n <- length(labels)
  boot_vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    boot_vals[b] <- metric_fn(labels[idx], scores[idx])
  }
  boot_vals <- boot_vals[!is.na(boot_vals)]
  quantile(boot_vals, c(0.025, 0.975))
}

cat("\n=== Final metrics with bootstrap 95% CIs (n_boot=2000) ===\n")
for (win in c("validate", "test")) {
  d <- fread(sprintf("phase8_out/phase8_%s_scored.csv", win))
  auroc <- calc_auc(d$became_dual_risk, d$predicted_risk)
  auprc <- calc_auprc(d$became_dual_risk, d$predicted_risk)
  auroc_ci <- boot_ci(d$became_dual_risk, d$predicted_risk, calc_auc)
  auprc_ci <- boot_ci(d$became_dual_risk, d$predicted_risk, calc_auprc)
  cat(sprintf("\n%s (n=%d):\n", win, sum(!is.na(d$predicted_risk))))
  cat(sprintf("  AUROC = %.3f [%.3f, %.3f]  (Phase 6 baseline: 0.590)\n", auroc, auroc_ci[1], auroc_ci[2]))
  cat(sprintf("  AUPRC = %.3f [%.3f, %.3f]  (Phase 6 baseline: 0.233, base rate: %.3f)\n",
              auprc, auprc_ci[1], auprc_ci[2], mean(d$became_dual_risk)))
}