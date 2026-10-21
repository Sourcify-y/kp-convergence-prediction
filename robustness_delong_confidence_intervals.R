suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pROC)
})

# AUROC + DeLong 95% CI for one scored table (rows with NA predicted_risk are dropped first)
score_file <- function(path, label) {
  if (!file.exists(path)) {
    cat(sprintf("[SKIP] %s — file not found: %s\n", label, path))
    return(NULL)
  }
  d <- read_csv(path, show_col_types = FALSE)

  n_total <- nrow(d)
  n_na    <- sum(is.na(d$predicted_risk))
  scored  <- d %>% filter(!is.na(predicted_risk))
  n_scored <- nrow(scored)
  n_pos <- sum(scored$became_dual_risk == 1, na.rm = TRUE)
  n_neg <- sum(scored$became_dual_risk == 0, na.rm = TRUE)

  if (n_pos == 0 || n_neg == 0) {
    cat(sprintf(
      "\n=== %s ===\nn_total=%d  NA=%d  scored=%d  pos=%d  neg=%d\n  Cannot compute AUROC/CI — need both classes present.\n",
      label, n_total, n_na, n_scored, n_pos, n_neg
    ))
    return(tibble(
      label = label, n_total = n_total, n_na = n_na, n_scored = n_scored,
      n_pos = n_pos, n_neg = n_neg, auroc = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_
    ))
  }

  roc_obj <- roc(scored$became_dual_risk, scored$predicted_risk, quiet = TRUE, direction = "<")
  auroc_val <- as.numeric(auc(roc_obj))
  ci_val <- ci.auc(roc_obj, method = "delong")

  cat(sprintf(
    "\n=== %s ===\nn_total=%d  NA=%d  scored=%d  pos=%d  neg=%d\n  AUROC = %.3f   95%% DeLong CI = [%.3f, %.3f]\n",
    label, n_total, n_na, n_scored, n_pos, n_neg, auroc_val, ci_val[1], ci_val[3]
  ))

  tibble(
    label = label, n_total = n_total, n_na = n_na, n_scored = n_scored,
    n_pos = n_pos, n_neg = n_neg, auroc = auroc_val,
    ci_lower = ci_val[1], ci_upper = ci_val[3]
  )
}

files <- list(
  list(path = "phase9_out/phase9_validate_holdout_scored.csv", label = "Lineage-holdout VALIDATE"),
  list(path = "phase9_out/phase9_test_holdout_scored.csv",     label = "Lineage-holdout TEST"),
  list(path = "phase9_out/phase9_validate_geo_scored.csv",     label = "Geo-holdout VALIDATE"),
  list(path = "phase9_out/phase9_test_geo_scored.csv",         label = "Geo-holdout TEST"),
  list(path = "phase8_out/phase8_sensitivity_check_pooled_validate.csv", label = "Sensitivity-check VALIDATE (3-seed pooled)"),
  list(path = "phase8_out/phase8_sensitivity_check_pooled_test.csv",     label = "Sensitivity-check TEST (3-seed pooled)")
)

results <- lapply(files, function(f) score_file(f$path, f$label))
results <- bind_rows(results)

cat("\n\n========== SUMMARY TABLE ==========\n")
print(results, n = Inf, width = Inf)

dir.create("phase9_out", showWarnings = FALSE)
write_csv(results, "phase9_out/delong_summary_all.csv")
cat("\nWritten to phase9_out/delong_summary_all.csv\n")