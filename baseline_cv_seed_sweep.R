suppressMessages(library(data.table))

dat <- fread("phase3_out/character_matrix.csv")
train_dat <- dat[collection_year <= 2017]

# most common country per ST, then bucket everything outside the top 5 as "Other"
st_info <- train_dat[, .(top_country = names(sort(table(country), decreasing = TRUE))[1]), by = base_ST]
top5 <- names(sort(table(st_info$top_country), decreasing = TRUE))[1:5]
st_info[, country_bucket := fifelse(top_country %in% top5, top_country, "Other")]

# rank-based AUROC (Mann-Whitney)
calc_auc <- function(y, score) {
  rk <- rank(score)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  (sum(rk[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# k-fold CV baseline: starting_state + country only
add_baseline <- function(task, k = 5) {
  d <- merge(task[, .(base_ST, starting_state, became_dual_risk, predicted_risk)],
             st_info[, .(base_ST, country_bucket)], by = "base_ST")
  d[, country_bucket := factor(country_bucket)]
  d[, starting_state := factor(starting_state)]
  fold_id <- sample(rep(1:k, length.out = nrow(d)))
  d[, baseline_pred := NA_real_]
  for (f in 1:k) {
    train_i <- which(fold_id != f)
    test_i <- which(fold_id == f)
    fit <- suppressWarnings(glm(became_dual_risk ~ starting_state + country_bucket,
                                data = d[train_i], family = binomial))
    d$baseline_pred[test_i] <- suppressWarnings(predict(fit, newdata = d[test_i], type = "response"))
  }
  d
}

# one-sided bootstrap p: how often model AUC <= baseline AUC
boot_pval <- function(y, model_score, base_score, n_boot = 1000) {
  n <- length(y)
  diffs <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    diffs[i] <- calc_auc(y[idx], model_score[idx]) - calc_auc(y[idx], base_score[idx])
  }
  mean(diffs <= 0, na.rm = TRUE)
}

for (split in c("validate", "test")) {
  scored <- fread(sprintf("phase8_out/phase8_%s_scored.csv", split))
  res <- data.table()
  for (s in 1:200) {
    set.seed(s)
    d <- add_baseline(scored)
    b_auc <- calc_auc(d$became_dual_risk, d$baseline_pred)
    m_auc <- calc_auc(d$became_dual_risk, d$predicted_risk)
    res <- rbind(res, data.table(seed = s, base_auc = b_auc, diff = m_auc - b_auc,
                                 p1 = boot_pval(d$became_dual_risk, d$predicted_risk, d$baseline_pred)))
  }
  # m_auc doesn't depend on the seed so the last one is fine to report
  cat(sprintf("\n%s (n=%d, model AUROC=%.3f)\n", toupper(split), nrow(d), m_auc))
  cat(sprintf("  baseline AUROC: mean=%.3f sd=%.3f min=%.3f max=%.3f\n",
              mean(res$base_auc), sd(res$base_auc), min(res$base_auc), max(res$base_auc)))
  cat(sprintf("  diff: mean=%.3f min=%.3f max=%.3f\n", mean(res$diff), min(res$diff), max(res$diff)))
  cat(sprintf("  one-sided p: median=%.3f | seeds with p<0.025: %.0f%% | seeds with p<0.00625: %.0f%%\n",
              median(res$p1), 100 * mean(res$p1 < 0.025), 100 * mean(res$p1 < 0.00625)))
}