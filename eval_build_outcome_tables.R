# eval_build_outcome_tables.R
options(warn = 1)
library(data.table)

dat <- fread("phase3_out/character_matrix.csv")
need <- c("base_ST", "collection_year",
          "AGly_acquired","Bla_acquired","Bla_ESBL_acquired","Bla_Carb_acquired",
          "Col_acquired","Fcyn_acquired","Flq_acquired","MLS_acquired","Phe_acquired",
          "Rif_acquired","Sul_acquired","Tet_acquired","Tgc_acquired","Tmt_acquired",
          "Yersiniabactin","Aerobactin","Salmochelin","Colibactin","rmpA","rmpA2")
missing_cols <- setdiff(need, names(dat))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse=", "),
                                   "\nActual columns: ", paste(names(dat), collapse=", "))

res_cols <- c("AGly_acquired","Bla_acquired","Bla_ESBL_acquired","Bla_Carb_acquired",
              "Col_acquired","Fcyn_acquired","Flq_acquired","MLS_acquired","Phe_acquired",
              "Rif_acquired","Sul_acquired","Tet_acquired","Tgc_acquired","Tmt_acquired")
vir_cols <- c("Yersiniabactin","Aerobactin","Salmochelin","Colibactin","rmpA","rmpA2")

# does each lineage have any resistance / virulence character in this window?
st_state <- function(win_dat) {
  win_dat[, .(
    has_resistance = as.integer(any(unlist(.SD[, res_cols, with = FALSE]) == 1, na.rm = TRUE)),
    has_virulence  = as.integer(any(unlist(.SD[, vir_cols, with = FALSE]) == 1, na.rm = TRUE))
  ), by = base_ST]
}

build_task <- function(train_end, win_start, win_end, out_path) {
  train_st <- st_state(dat[collection_year <= train_end])
  win_dat <- dat[collection_year >= win_start & collection_year <= win_end]
  win_st <- st_state(win_dat)
  setnames(win_st, c("has_resistance","has_virulence"), c("has_resistance_W","has_virulence_W"))

  # inner join: lineages with no genomes in the window drop out (same as Phase 6)
  task <- merge(train_st, win_st, by = "base_ST")
  task <- task[(has_resistance == 1) != (has_virulence == 1)]
  task[, became_dual_risk := as.integer(has_resistance_W == 1 & has_virulence_W == 1)]
  task[, starting_state := fifelse(has_resistance == 1, "resistance_only", "virulence_only")]

  out <- task[, .(base_ST, starting_state, became_dual_risk)]
  fwrite(out, out_path)
  cat(sprintf("%s: %d single-risk lineages, %.1f%% became dual-risk\n",
              out_path, nrow(out), 100 * mean(out$became_dual_risk)))
}

build_task(2017, 2018, 2020, "phase8_out/phase8_validate_task_table.csv")  # comparable to Phase 6's 0.590/0.233
build_task(2017, 2021, 2022, "phase8_out/phase8_test_task_table.csv")