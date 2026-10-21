#!/usr/bin/env Rscript
#
# Phase 4: tree building + ancestral state reconstruction
# Kp hv/CRKP convergence project
#
# Builds a phylogenetic tree for the 8,318-genome Mash-deduplicated, usable dataset, then
# reconstructs ancestral states for the 20 locked characters (14 resistance + 6 virulence) so
# Phase 4's acquisition-event definition (character absent at the ancestral node, present at
# the descendant) can be evaluated.
#
# Design choices:
#   1. The tree is built by neighbor-joining on the Mash distance matrix already computed for
#      the dedup step, not from a fresh core-genome alignment. This reuses the dedup distances
#      and avoids a multi-day alignment across 8,318 genomes.
#   2. ASR uses PARSIMONY (Fitch, via phangorn), not maximum likelihood. Mash-derived branch
#      lengths are a distance proxy, not substitutions/site, so an ML model that depends on
#      calibrated branch lengths (e.g. ace() with ER/SYM) is on shakier ground here. Parsimony
#      also gives a natural per-node ambiguity signal, useful for the confidence-threshold
#      question in Phase 4 Section 3.1. ML ASR is a documented sensitivity-check TODO (Phase
#      9-style) and is not run by default.
#   3. NJ trees are unrooted, so the tree is midpoint-rooted, because the directional
#      "ancestral absent -> descendant present" logic needs a root. This is a real assumption
#      to flag at sign-off: midpoint rooting can be wrong when substitution rates are very
#      uneven across lineages, which is plausible given Kp's recombination rate (Phase 4 open
#      item #1). Treat it as provisional, not settled.
#
# Requires (CRAN): ape, phangorn, optparse, data.table
#   install.packages(c("ape","phangorn","optparse","data.table"))
#
# Usage:
#   Rscript tree_build_and_acquisition_events.R \
#     --distances kp_pull_20k/mash_test/distances_full.tsv \
#     --representatives kp_pull_20k/mash_test/representative_genomes.txt \
#     --manifest kp_pull_20k/03_download_manifest.csv \
#     --matrix phase3_out/character_matrix.csv \
#     --outdir phase4_out

suppressMessages({
  library(optparse)
  library(data.table)
  library(ape)
  library(phangorn)
})

# ---------------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--distances", type = "character",
              help = "Path to mash triangle output (lower-triangular, full 15,764-genome matrix)"),
  make_option("--representatives", type = "character",
              help = "Path to representative_genomes.txt (9,108 paths, Mash-dedup output)"),
  make_option("--manifest", type = "character",
              help = "Path to 03_download_manifest.csv mapping genome file paths -> assembly_accession"),
  make_option("--matrix", type = "character",
              help = "Path to phase3_out/character_matrix.csv (8,318 x 24: 20 characters + genome_id, collection_year, base_ST, country)"),
  make_option("--outdir", type = "character", default = "phase4_out",
              help = "Output directory [default %default]"),
  make_option("--tree", type = "character", default = NULL,
              help = "Optional: path to a previously saved tree (e.g. phase4_out/kp_dedup_nj_tree.nwk). Skips the NJ build (~6-7 min at 8,318 genomes) and goes straight to ASR. --distances/--representatives are still used to load and align the character matrix's genome set."),
  make_option("--min_confidence", type = "double", default = 0.9,
              help = "Minimum parsimony-state proportion at an internal node to call it non-ambiguous [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("distances", "representatives", "manifest", "matrix")) {
  if (is.null(opt[[req]])) stop(sprintf("--%s is required", req))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# ---------------------------------------------------------------------------
# Step 1: load the character matrix (defines the exact genome set + traits)
# ---------------------------------------------------------------------------
log_msg("Loading character matrix from %s", opt$matrix)
char_mat <- fread(opt$matrix)

id_cols <- intersect(c("genome_id", "assembly_accession", "collection_year",
                       "base_ST", "country"), names(char_mat))
if (!"genome_id" %in% names(char_mat) && "assembly_accession" %in% names(char_mat)) {
  setnames(char_mat, "assembly_accession", "genome_id")
}
if (!"genome_id" %in% names(char_mat)) {
  stop("character_matrix.csv must have a genome_id (or assembly_accession) column to match against tree tips")
}

char_cols <- setdiff(names(char_mat), id_cols)
log_msg("Detected %d character columns (expect 20): %s",
        length(char_cols), paste(char_cols, collapse = ", "))
if (length(char_cols) != 20) {
  warning(sprintf(
    "Expected exactly 20 character columns (14 resistance + 6 virulence), found %d. Check column detection before trusting downstream output.",
    length(char_cols)))
}

# ---------------------------------------------------------------------------
# Step 2: genome_id -> tree-tip-label mapping via the manifest
# ---------------------------------------------------------------------------
log_msg("Loading manifest from %s", opt$manifest)
manifest <- fread(opt$manifest)
# the manifest's column names aren't pinned anywhere, so look for them by pattern
path_col <- names(manifest)[grepl("path|file", names(manifest), ignore.case = TRUE)][1]
acc_col  <- names(manifest)[grepl("accession", names(manifest), ignore.case = TRUE)][1]
if (is.na(path_col) || is.na(acc_col)) {
  stop("Could not identify path/accession columns in manifest — inspect 03_download_manifest.csv headers and adjust this script's column detection.")
}
manifest[[path_col]] <- gsub("\r", "", manifest[[path_col]])  # guard against the CRLF bug found earlier in this project

log_msg("Loading representative genome list from %s", opt$representatives)
reps <- fread(opt$representatives, header = FALSE)[[1]]
reps <- gsub("\r", "", reps)

rep_acc <- manifest[[acc_col]][match(reps, manifest[[path_col]])]
if (any(is.na(rep_acc))) {
  n_missing <- sum(is.na(rep_acc))
  warning(sprintf("%d/%d representative genome paths did not match the manifest — check for path formatting drift before trusting the tip set.",
                  n_missing, length(reps)))
  rep_acc <- rep_acc[!is.na(rep_acc)]
}

# Keep genomes present in BOTH the Mash-dedup representative set and the usable character
# matrix. The matrix is already filtered to usable + dedup genomes (Phase 4 Sec 4), so this
# should be close to a no-op. If it isn't, investigate rather than silently discarding.
keep_ids <- intersect(rep_acc, char_mat$genome_id)
log_msg("%d genomes in both the Mash-dedup representative set and the character matrix (expect close to 8,318)",
        length(keep_ids))
if (length(keep_ids) < 0.95 * nrow(char_mat)) {
  warning("Fewer than 95% of character_matrix.csv genomes matched to the representative set — investigate before proceeding to ASR.")
}

char_mat <- char_mat[genome_id %in% keep_ids]

# ---------------------------------------------------------------------------
# Step 3: load the full Mash distance matrix, subset to the kept genomes
# ---------------------------------------------------------------------------
log_msg("Loading Mash distance matrix from %s (this is a large file, ~1.26GB per earlier run — allow time/RAM)", opt$distances)
# mash triangle output: first line = N, then a lower-triangular matrix with genome names as
# row labels and no header row, so it's read by hand instead of assuming a clean CSV
con <- file(opt$distances, "r")
n_genomes <- as.integer(readLines(con, n = 1))
names_vec <- character(n_genomes)
dist_list <- vector("list", n_genomes)
for (i in seq_len(n_genomes)) {
  line <- readLines(con, n = 1)
  parts <- strsplit(line, "\t")[[1]]
  names_vec[i] <- parts[1]
  if (length(parts) > 1) {
    dist_list[[i]] <- as.numeric(parts[-1])
  } else {
    dist_list[[i]] <- numeric(0)
  }
}
close(con)
names_vec <- gsub("\r", "", names_vec)

# rebuild a full symmetric matrix
full_mat <- matrix(0, n_genomes, n_genomes, dimnames = list(names_vec, names_vec))
for (i in seq_len(n_genomes)) {
  if (length(dist_list[[i]]) > 0) {
    full_mat[i, seq_len(length(dist_list[[i]]))] <- dist_list[[i]]
    full_mat[seq_len(length(dist_list[[i]])), i] <- dist_list[[i]]
  }
}

# the names in the distance matrix are file paths, so map them back to accessions
mat_acc <- manifest[[acc_col]][match(names_vec, manifest[[path_col]])]
rownames(full_mat) <- mat_acc
colnames(full_mat) <- mat_acc

sub_ids <- intersect(keep_ids, mat_acc)
log_msg("%d genomes matched between distance matrix and kept genome set — subsetting distance matrix", length(sub_ids))
sub_mat <- full_mat[sub_ids, sub_ids]
char_mat <- char_mat[genome_id %in% sub_ids]
char_mat <- char_mat[match(sub_ids, genome_id)]  # same row order as the distance matrix

# ---------------------------------------------------------------------------
# Step 4: build the NJ tree and midpoint-root it
# ---------------------------------------------------------------------------
log_msg("Building neighbor-joining tree on %d genomes", length(sub_ids))
if (!is.null(opt$tree)) {
  log_msg("--tree supplied (%s) — loading saved tree instead of rebuilding (skips the ~6-7 min NJ step)", opt$tree)
  tree <- read.tree(opt$tree)
  if (!setequal(tree$tip.label, sub_ids)) {
    warning(sprintf(
      "Loaded tree's tip set does not exactly match the current genome set (%d tips loaded vs %d expected) — this can happen if --matrix or --representatives changed since the tree was built. Proceeding with the loaded tree's own tip set, but double check this is intentional.",
      length(tree$tip.label), length(sub_ids)))
    char_mat <- char_mat[genome_id %in% tree$tip.label]
    char_mat <- char_mat[match(tree$tip.label, genome_id)]
  }
} else {
  tree <- nj(as.dist(sub_mat))
  log_msg("Midpoint-rooting tree (NJ trees are unrooted; ASR direction requires a root — see header notes on this assumption)")
  tree <- midpoint(tree)
  tree <- ladderize(tree)

  tree_path <- file.path(opt$outdir, "kp_dedup_nj_tree.nwk")
  write.tree(tree, tree_path)
  log_msg("Wrote tree to %s", tree_path)
}

# ---------------------------------------------------------------------------
# Step 5: parsimony ASR per character
# ---------------------------------------------------------------------------
log_msg("Running Fitch-parsimony ASR for %d characters", length(char_cols))

acquisition_events <- list()
asr_by_char <- list()

for (ch in char_cols) {
  log_msg("  -> %s", ch)
  states <- as.character(char_mat[[ch]])
  names(states) <- char_mat$genome_id

  # phyDat object for a binary (0/1) character
  phy <- phyDat(matrix(states[tree$tip.label], ncol = 1,
                       dimnames = list(tree$tip.label, ch)),
                type = "USER", levels = c("0", "1"))

  anc <- tryCatch(
    ancestral.pars(tree, phy, type = "MPR"),
    error = function(e) {
      warning(sprintf("ASR failed for character %s: %s", ch, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(anc)) next

  # anc is a list over all nodes (tips + internal), each a probability vector over the phyDat
  # levels. For binary parsimony it's typically 0, 1 or 0.5 (ambiguous). phangorn doesn't
  # attach dimnames to these matrices; columns follow the "levels" order passed to phyDat()
  # above (c("0", "1")), so the "1" column has to be found by position, not by name
  # (indexing by name fails with "no 'dimnames' attribute for array").
  lvls <- attr(phy, "levels")
  state1_col <- which(lvls == "1")
  if (length(state1_col) != 1) {
    warning(sprintf("Could not uniquely locate the '1' state column for character %s (levels: %s) — skipping",
                    ch, paste(lvls, collapse = ",")))
    next
  }
  node_labels <- names(anc)
  n_nodes_total <- length(anc)
  expected_nodes <- length(tree$tip.label) + tree$Nnode
  if (is.null(node_labels)) {
    log_msg("    (ancestral.pars() output has no names — indexing by node-number position directly; %d nodes in list, %d expected from tree)",
            n_nodes_total, expected_nodes)
  }
  prob_1 <- vapply(anc, function(x) x[1, state1_col], numeric(1))
  is_confident <- prob_1 >= opt$min_confidence | prob_1 <= (1 - opt$min_confidence)
  # Don't break exact ties (prob_1 == 0.5) toward "1". A real Fitch tie means the ancestral
  # state is genuinely unresolved at that node. Defaulting ties to "1" would structurally
  # exclude every ambiguous node from being an event's ancestral node (an event needs
  # anc_state == "0"), which hides real ambiguity instead of reporting it. This was caught
  # because the confidence diagnostic read as a suspicious, uninformative 100% across all
  # recorded events.
  called_state <- ifelse(prob_1 > 0.5, "1",
                         ifelse(prob_1 < 0.5, "0", NA_character_))
  n_ambiguous <- sum(is.na(called_state))
  log_msg("    %d/%d nodes genuinely ambiguous (exact parsimony ties) for %s",
          n_ambiguous, length(called_state), ch)

  names(called_state) <- as.character(seq_along(called_state))
  asr_by_char[[ch]] <- called_state
  # tree$edge: ancestor node -> descendant node for each edge, in ape numbering
  # (tips 1..n, internal n+1..)
  n_tip <- length(tree$tip.label)
  for (e in seq_len(nrow(tree$edge))) {
    anc_node <- tree$edge[e, 1]
    desc_node <- tree$edge[e, 2]
    # ancestral.pars()'s list has NO names: position i is node number i (same convention
    # plotAnc() uses via anc[[node]]), so index by node number, not by a name string
    if (anc_node > n_nodes_total || desc_node > n_nodes_total) next

    anc_state <- called_state[anc_node]
    desc_state <- called_state[desc_node]
    anc_conf <- is_confident[anc_node]

    # Skip edges where either end is unresolved (exact parsimony tie). That's real ambiguity,
    # so it isn't treated as a non-event.
    if (is.na(anc_state) || is.na(desc_state)) next

    # candidate acquisition: absent -> present on this edge
    if (length(anc_state) && length(desc_state) &&
        anc_state == "0" && desc_state == "1") {

      desc_is_tip <- desc_node <= n_tip
      acquisition_events[[length(acquisition_events) + 1]] <- data.table(
        character = ch,
        ancestral_node = anc_node,
        descendant_node = desc_node,
        descendant_is_tip = desc_is_tip,
        descendant_genome_id = if (desc_is_tip) tree$tip.label[desc_node] else NA_character_,
        ancestral_state_confident = anc_conf,
        edge_length = tree$edge.length[e]
      )
    }
  }
}

asr_path <- file.path(opt$outdir, "asr_states_by_character.rds")
saveRDS(asr_by_char, asr_path)
log_msg("Wrote per-character ASR state (named vector per character, node-id keyed) to %s", asr_path)

events_dt <- rbindlist(acquisition_events, fill = TRUE)
if (nrow(events_dt) == 0) {
  stop(paste0(
    "Zero candidate acquisition events detected across all ", length(char_cols),
    " characters. Given how common these traits are in your data, this almost ",
    "certainly means a bug in the node-matching logic above, not a genuine biological ",
    "result — do not proceed to Section 3 exclusion filtering on this output."
  ))
}

# add collection_year for tip-level events, for the temporal filtering later (Phase 4
# Section 3.1: the ancestral-presence / sampling-artifact check benefits from knowing how
# early or late a flagged tip was collected)
if ("collection_year" %in% names(char_mat)) {
  events_dt[descendant_is_tip == TRUE,
            collection_year := char_mat$collection_year[match(descendant_genome_id, char_mat$genome_id)]]
}

events_path <- file.path(opt$outdir, "candidate_acquisition_events.csv")
fwrite(events_dt, events_path)
log_msg("Wrote %d candidate acquisition-events (pre-exclusion-filter) to %s", nrow(events_dt), events_path)

# ---------------------------------------------------------------------------
# Step 6: summary (not a substitute for the Section 3 exclusion logic)
# ---------------------------------------------------------------------------
log_msg(paste0(
  "Done. IMPORTANT: this script only implements the raw ancestral-absent -> ",
  "descendant-present detection (Phase 4 Sec 2). It does NOT yet apply ",
  "Section 3's four exclusions (ancestral-presence/sampling artifact, ",
  "loss/reappearance, HGT/recombination tagging, assembly artifact) — ",
  "those need dedicated logic layered on top of candidate_acquisition_events.csv ",
  "before any of this is used for Phase 6-9 counts."
))