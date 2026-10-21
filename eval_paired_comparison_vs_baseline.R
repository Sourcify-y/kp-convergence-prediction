# eval_paired_comparison_vs_baseline.R
library(data.table)
set.seed(42)

res_chars <- c("AGly_acquired","Bla_acquired","Bla_ESBL_acquired","Bla_Carb_acquired",
               "Col_acquired","Fcyn_acquired","Flq_acquired","MLS_acquired","Phe_acquired",
               "Rif_acquired","Sul_acquired","Tet_acquired","Tgc_acquired","Tmt_acquired")
vir_chars <- c("Yersiniabactin","Aerobactin","Salmochelin","Colibactin","rmpA","rmpA2")

dat <- fread("phase3_out/character_matrix.csv")
if (!"country" %in% names(dat)) stop("No 'country' column. Actual columns: ", paste(names(dat), collapse=", "))
train_dat <- dat[collection_year <= 2017]

st_info <- train_dat[, .(top_country = names(sort(table(country), decreasing = TRUE))[1]),
                     by = base_ST]
# bucket rare countries so CV folds don't hit unseen levels
top5 <- names(sort(table(st_info$top_country), decreasing = TRUE))[1:5]
st_info[, country_bucket := fifelse(top_country %in% top5, top_country, "Other")]

add_baseline <- function(task, k = 5) {
  d <- merge(task[, .(base_ST, starting_state, became_dual_risk)],
             st_info[, .(base_ST, country_bucket)], by = "base_ST")
  d[, country_bucket := factor(country_bucket)]
  d[, starting_state := factor(starting_state)]
  n <- nrow(d)
  fold_id <- sample(rep(1:k, length.out = n))
  d[, baseline_pred := NA_real_]
  for (f in 1:k) {
    train_i <- which(fold_id != f); test_i <- which(fold_id == f)
    fit <- glm(became_dual_risk ~ starting_state + country_bucket, data = d[train_i], family = binomial)
    d$baseline_pred[test_i] <- predict(fit, newdata = d[test_i], type = "response")
  }
  d
}

val_scored <- fread("phase8_out/phase8_validate_scored.csv")
test_scored <- fread("phase8_out/phase8_test_scored.csv")

boot_auc_diff <- function(labels, score_a, score_b, n_boot = 2000) {
  n <- length(labels)
  diffs <- numeric(n_boot)
  calc_auc <- function(l, s) { r <- rank(s); n1 <- sum(l==1); n0 <- sum(l==0)
    (sum(r[l==1]) - n1*(n1+1)/2) / (n1*n0) }
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    diffs[b] <- calc_auc(labels[idx], score_a[idx]) - calc_auc(labels[idx], score_b[idx])
  }
  list(diff = mean(diffs), ci = quantile(diffs, c(0.025, 0.975)),
       p_one_sided = mean(diffs <= 0))
}

for (win in c("validate", "test")) {
  scored <- if (win == "validate") val_scored else test_scored
  base_df <- add_baseline(scored)
  mrg <- merge(scored[, .(base_ST, became_dual_risk, predicted_risk)],
               base_df[, .(base_ST, baseline_pred)], by = "base_ST")
  out <- boot_auc_diff(mrg$became_dual_risk, mrg$predicted_risk, mrg$baseline_pred)
  cat(sprintf("\n%s (n=%d): HyperTraPS - baseline AUROC diff = %.3f, 95%% CI [%.3f, %.3f], P(HyperTraPS <= baseline) = %.3f\n",
              win, nrow(mrg), out$diff, out$ci[1], out$ci[2], out$p_one_sided))
}

# AUPRC (step-wise average precision) and a generic paired bootstrap for any metric
calc_auprc <- function(l, s) {
  o <- order(s, decreasing = TRUE)
  l <- l[o]
  tp <- cumsum(l); fp <- cumsum(1 - l)
  precision <- tp / (tp + fp); recall <- tp / sum(l)
  sum(diff(c(0, recall)) * precision)
}

boot_metric_diff <- function(labels, score_a, score_b, metric_fn, n_boot = 10000, seed = 42) {
  set.seed(seed)
  n <- length(labels)
  diffs <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    diffs[b] <- metric_fn(labels[idx], score_a[idx]) - metric_fn(labels[idx], score_b[idx])
  }
  list(diff = mean(diffs), ci = quantile(diffs, c(0.025, 0.975)), p_one_sided = mean(diffs <= 0))
}

calc_auc <- function(l, s) { r <- rank(s); n1 <- sum(l==1); n0 <- sum(l==0)
  (sum(r[l==1]) - n1*(n1+1)/2) / (n1*n0) }

tab <- data.table()
for (win in c("validate", "test")) {
  scored <- if (win == "validate") val_scored else test_scored
  base_df <- add_baseline(scored)
  mrg <- merge(scored[, .(base_ST, became_dual_risk, predicted_risk)],
               base_df[, .(base_ST, baseline_pred)], by = "base_ST")

  for (metric_name in c("AUROC", "AUPRC")) {
    metric_fn <- if (metric_name == "AUROC") calc_auc else calc_auprc
    out <- boot_metric_diff(mrg$became_dual_risk, mrg$predicted_risk, mrg$baseline_pred, metric_fn)
    tab <- rbind(tab, data.table(
      window = win, metric = metric_name, n = nrow(mrg),
      diff = round(out$diff, 3),
      ci_low = round(out$ci[1], 3), ci_high = round(out$ci[2], 3),
      p_no_improvement = round(out$p_one_sided, 3)
    ))
  }
}

cat("\n=== Full paired comparison table (n_boot=10000, seed=42) — all results, none excluded ===\n")
print(tab)
fwrite(tab, "phase8_out/phase8_paired_comparison_full.csv")

brier_score <- function(labels, scores) mean((labels - scores)^2, na.rm = TRUE)

for (win in c("validate", "test")) {
  d <- fread(sprintf("phase8_out/phase8_%s_scored.csv", win))
  cat(sprintf("%s Brier score = %.3f (n=%d)\n", win, brier_score(d$became_dual_risk, d$predicted_risk), sum(!is.na(d$predicted_risk))))
}

recal_brier <- function(labels, scores) {
  d <- data.frame(labels = labels, scores = scores)
  fit <- glm(labels ~ scores, data = d, family = binomial)
  recal <- predict(fit, type = "response")
  mean((labels - recal)^2, na.rm = TRUE)
}

for (win in c("validate", "test")) {
  d <- fread(sprintf("phase8_out/phase8_%s_scored.csv", win))
  cat(sprintf("%s recalibrated Brier = %.3f\n", win, recal_brier(d$became_dual_risk, d$predicted_risk)))
}