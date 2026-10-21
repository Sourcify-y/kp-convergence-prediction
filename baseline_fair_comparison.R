library(data.table)

calc_auc <- function(y, score) {
  rk <- rank(score)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  (sum(rk[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

dat <- fread("phase3_out/character_matrix.csv")
train_dat <- dat[collection_year <= 2017]
st_info <- train_dat[, .(top_country = names(sort(table(country), decreasing = TRUE))[1]), by = base_ST]
top5 <- names(sort(table(st_info$top_country), decreasing = TRUE))[1:5]
st_info[, country_bucket := fifelse(top_country %in% top5, top_country, "Other")]

load_task <- function(split) {
  scored <- fread(sprintf("phase8_out/phase8_%s_scored.csv", split))
  merge(scored[, .(base_ST, starting_state, became_dual_risk, predicted_risk)],
        st_info[, .(base_ST, country_bucket)], by = "base_ST")
}
val_df <- load_task("validate")
test_df <- load_task("test")

# shared factor levels so val and test line up
cb_levels <- sort(unique(c(val_df$country_bucket, test_df$country_bucket)))
ss_levels <- sort(unique(c(val_df$starting_state, test_df$starting_state)))

set_levels <- function(d) {
  d <- copy(d)
  d[, country_bucket := factor(country_bucket, levels = cb_levels)]
  d[, starting_state := factor(starting_state, levels = ss_levels)]
  d
}
val_df <- set_levels(val_df)
test_df <- set_levels(test_df)

fit_predict <- function(train_d, test_d) {
  fit <- suppressWarnings(glm(became_dual_risk ~ starting_state + country_bucket,
                              data = train_d, family = binomial))
  suppressWarnings(as.numeric(predict(fit, newdata = test_d, type = "response")))
}

report_auc <- function(label, d, base) {
  keep_i <- !is.na(base)
  d <- d[keep_i]
  base <- base[keep_i]
  set.seed(42)
  n <- nrow(d)
  diffs <- replicate(10000, {
    i <- sample(n, n, replace = TRUE)
    calc_auc(d$became_dual_risk[i], d$predicted_risk[i]) - calc_auc(d$became_dual_risk[i], base[i])
  })
  cat(sprintf("%s: n=%d model=%.3f baseline=%.3f | diff=%.3f, 95%% CI [%.3f, %.3f], p(diff<=0)=%.3f\n",
      label, n, calc_auc(d$became_dual_risk, d$predicted_risk), calc_auc(d$became_dual_risk, base),
      mean(diffs), quantile(diffs, 0.025), quantile(diffs, 0.975), mean(diffs <= 0)))
}

cat("\n--- In-sample baseline (optimistic for baseline, conservative for model) ---\n")
report_auc("VALIDATE in-sample", val_df, fit_predict(val_df, val_df))
report_auc("TEST in-sample    ", test_df, fit_predict(test_df, test_df))
cat("\n--- Out-of-window baseline (fit on VALIDATE, score TEST) ---\n")
report_auc("TEST fit-on-VALIDATE", test_df, fit_predict(val_df, test_df))