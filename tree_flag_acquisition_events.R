# ============================================================
# Phase 4 section 3: acquisition-event exclusion flags
#   Exclusion 2 (loss/reappearance, singleton/cluster framing)
#   -> produces candidate_acquisition_events_flagged_v2.csv
# ============================================================
#
# PROVENANCE NOTE: the console session that produced the original v2 table was not saved as
# a script. This file is a reconstruction. The reconstructed pieces were checked against the
# shipped v2 table:
#   - subtree_size = number of nodes within max_branch_dist = 0.0024 (cumulative branch
#     length, negative edges clamped to 0) of the event's descendant node: matches v2 for
#     7,393 / 7,393 events
#   - is_singleton = (subtree_size == 1): matches 7,393 / 7,393 (5,292 singletons, 2,101
#     cluster members)
#   - subtree_reverts is NA exactly for singletons: 7,393 / 7,393
# is_reacquisition and subtree_reverts were confirmed by running this script on the full
# event table: 0 mismatches against the shipped v2 table in all four derived columns
# (subtree_size, is_singleton, is_reacquisition, subtree_reverts).
#
# Definitions
#   is_reacquisition : on the root-to-ancestral-node path, the trait goes 1 -> 0 at least
#                      once (present upstream, lost, then regained). Computed for all events.
#   is_singleton     : no other node within max_branch_dist below the event's descendant node
#                      (no downstream shared ancestry, so no independent evidence).
#   subtree_reverts  : among non-singleton events only, whether any node within
#                      max_branch_dist is in state 0. NA for singletons (undetermined, not FALSE).
# Rationale: Boyko et al. (singleton vs cluster events) and Gascuel & Steel 2020 (ASR power
# decays away from well-sampled tips). max_branch_dist = 0.0024 was chosen at the plateau of a
# calibration sweep (checkable fraction ~28.4% vs a ~29.1% ceiling at 0.0096).
#
# Inputs (produced by tree_build_and_acquisition_events.R):
#   phase4_out/candidate_acquisition_events.csv
#   phase4_out/kp_dedup_nj_tree.nwk
#   phase4_out/asr_states_by_character.rds
# ============================================================

library(ape)
library(phangorn)

# ---------------------------------------------------------------
# tree-navigation helpers
# ---------------------------------------------------------------

get_parent_node <- function(tree, node) {
  edge <- tree$edge
  hit <- which(edge[, 2] == node)
  if (length(hit) == 0) return(NA_integer_)  # node is the root
  edge[hit, 1]
}

get_children <- function(tree, node) {
  tree$edge[tree$edge[, 1] == node, 2]
}

get_sibling_node <- function(tree, parent, child) {
  kids <- get_children(tree, parent)
  sib <- kids[kids != child]
  if (length(sib) != 1) {
    stop(sprintf(
      paste0("Expected exactly one sibling of node %d under parent %d, ",
             "found %d — tree may not be strictly bifurcating at this node."),
      child, parent, length(sib)
    ))
  }
  sib
}

count_tips_in_clade <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) return(1L)  # node itself is a tip
  desc <- Descendants(tree, node, type = "tips")[[1]]
  length(desc)
}

# Minimum cumulative branch length from `node` to any tip in its own subtree. A distance-based
# (not topology-count) proxy for how well evidenced this part of the tree is, robust to uneven
# sampling depth across lineages (unlike raw sibling tip count).
min_branch_dist_to_tip <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) return(0)  # node is itself a tip

  dists <- setNames(NA_real_, as.character(node))
  dists[as.character(node)] <- 0
  frontier <- node
  min_tip_dist <- Inf

  while (length(frontier) > 0) {
    next_frontier <- integer(0)
    for (f in frontier) {
      kids <- get_children(tree, f)
      for (k in kids) {
        el <- tree$edge.length[tree$edge[, 1] == f & tree$edge[, 2] == k]
        el <- max(el, 0)  # clamp negative NJ branch lengths, same as elsewhere
        d <- dists[as.character(f)] + el
        dists[as.character(k)] <- d
        if (k <= n_tip) {
          min_tip_dist <- min(min_tip_dist, d)
        } else {
          next_frontier <- c(next_frontier, k)
        }
      }
    }
    frontier <- next_frontier
  }
  min_tip_dist
}

root_to_node_path <- function(tree, node) {
  # node IDs from the root down to `node`, both ends included, in order
  path <- node
  cur <- node
  repeat {
    p <- get_parent_node(tree, cur)
    if (is.na(p)) break
    path <- c(p, path)
    cur <- p
  }
  path
}

get_subtree_all_nodes <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) return(node)
  c(node, Descendants(tree, node, type = "all")[[1]])
}

# Nodes within `max_branch_dist` cumulative branch length of `node`, restricted to the node's
# own descendant clade (BFS over ape edge lengths).
#
# Negative-edge note: this NJ tree has 37/16634 (~0.22%) negative branch lengths (a known NJ
# artifact, not a pipeline bug), up to ~4.1e-4 in magnitude, the same order as the smallest
# max_branch_dist candidates. Left as is, a negative edge breaks the assumption that cumulative
# distance only grows moving away from `node`, which would let a node farther down the tree
# look closer than a nearer one and silently corrupt which nodes count as in range. So they're
# clamped to 0 here and a one-time warning logs the count. This is a deliberate, documented
# methodological choice, not a silent patch.
get_subtree_within_distance <- function(tree, node, max_branch_dist) {
  if (is.null(max_branch_dist)) {
    stop("max_branch_dist is unset — calibrate against your tree's branch-length distribution before running exclusion 2 (see header notes).")
  }

  n_negative <- sum(tree$edge.length < 0)
  if (n_negative > 0 && !isTRUE(getOption("phase4_neg_edge_warned"))) {
    warning(sprintf(
      "%d/%d edges (%.2f%%) have negative branch length (NJ artifact) — clamped to 0 for the exclusion-2 distance-limited BFS. Log this in spec.md as a stated assumption.",
      n_negative, length(tree$edge.length), 100 * n_negative / length(tree$edge.length)
    ))
    options(phase4_neg_edge_warned = TRUE)  # once per session, not once per event
  }

  all_desc <- get_subtree_all_nodes(tree, node)
  if (length(all_desc) == 1) return(all_desc)  # a tip, no descendants

  dists <- setNames(rep(NA_real_, length(all_desc)), as.character(all_desc))
  dists[as.character(node)] <- 0
  frontier <- node
  while (length(frontier) > 0) {
    next_frontier <- integer(0)
    for (f in frontier) {
      kids <- get_children(tree, f)
      for (k in kids) {
        el <- tree$edge.length[tree$edge[, 1] == f & tree$edge[, 2] == k]
        el <- max(el, 0)  # clamp negative NJ branch lengths, see note above
        d <- dists[as.character(f)] + el
        dists[as.character(k)] <- d
        if (d <= max_branch_dist) next_frontier <- c(next_frontier, k)
      }
    }
    frontier <- next_frontier
  }
  as.integer(names(dists)[!is.na(dists) & dists <= max_branch_dist])
}

# ---------------------------------------------------------------
# Exclusion 1: ancestral-presence / sampling artifact
# ---------------------------------------------------------------

flag_sampling_artifact <- function(event, tree, asr_states, min_tip_dist_threshold = 0.00015) {
  parent <- event$ancestral_node
  child  <- event$descendant_node
  state_vec <- asr_states[[event$character]]

  grandparent <- get_parent_node(tree, parent)
  grandparent_state <- NA_integer_
  if (!is.na(grandparent)) {
    grandparent_state <- state_vec[as.character(grandparent)]
  }

  sibling <- get_sibling_node(tree, parent, child)
  sibling_tip_count <- count_tips_in_clade(tree, sibling)
  sibling_min_tip_dist <- min_branch_dist_to_tip(tree, sibling)

  list(
    sibling_node = sibling,
    sibling_tip_count = sibling_tip_count,
    sibling_min_tip_dist = sibling_min_tip_dist,
    low_confidence_ancestral = sibling_min_tip_dist > min_tip_dist_threshold,
    grandparent_node = grandparent,
    grandparent_state = grandparent_state,
    possible_exclusion2_overlap = !is.na(grandparent_state) && grandparent_state == 1
  )
}

# ---------------------------------------------------------------
# Exclusion 2: loss / reappearance
# ---------------------------------------------------------------

count_state_changes <- function(state_vec, path_nodes) {
  states <- state_vec[as.character(path_nodes)]
  states <- states[!is.na(states)]
  if (length(states) < 2) {
    return(list(total = 0L, upstream_1_to_0_count = 0L, upstream_0_to_1_count = 0L))
  }
  diffs <- diff(states)
  list(
    total = sum(diffs != 0),
    upstream_1_to_0_count = sum(diffs == -1),
    upstream_0_to_1_count = sum(diffs == 1)
  )
}

any_reversion_in_subtree <- function(state_vec, subtree_nodes) {
  states <- state_vec[as.character(subtree_nodes)]
  states <- states[!is.na(states)]
  any(states == 0)
}

flag_loss_reappearance <- function(event, tree, asr_states, max_branch_dist = NULL) {
  state_vec <- asr_states[[event$character]]

  # full root-to-parent path (strictly upstream of the acquisition itself)
  full_path <- root_to_node_path(tree, event$ancestral_node)
  transitions <- count_state_changes(state_vec, full_path)

  subtree_nodes <- get_subtree_within_distance(tree, event$descendant_node, max_branch_dist)
  reverts <- any_reversion_in_subtree(state_vec, subtree_nodes)

  list(
    is_reacquisition = transitions$upstream_1_to_0_count > 0,
    upstream_transition_count = transitions$total,
    subtree_reverts = reverts
  )
}


# ---------------------------------------------------------------
# driver: writes the v2 table
# ---------------------------------------------------------------
V2_COLUMNS <- c("character", "ancestral_node", "descendant_node",
                "descendant_is_tip", "descendant_genome_id",
                "ancestral_state_confident", "edge_length",
                "collection_year", "subtree_size", "is_singleton",
                "is_reacquisition", "subtree_reverts")

run_final_exclusions <- function(
    events_csv_path     = "phase4_out/candidate_acquisition_events.csv",
    tree_path           = "phase4_out/kp_dedup_nj_tree.nwk",
    asr_states_rds_path = "phase4_out/asr_states_by_character.rds",
    max_branch_dist     = 0.0024,
    out_csv_path        = "phase4_out/v2_rebuilt.csv") {

  events <- read.csv(events_csv_path, stringsAsFactors = FALSE)
  need <- setdiff(V2_COLUMNS, c("subtree_size", "is_singleton",
                                "is_reacquisition", "subtree_reverts"))
  if (!all(need %in% names(events))) {
    stop("events file is missing columns: ",
         paste(setdiff(need, names(events)), collapse = ", "),
         "\nActual columns: ", paste(names(events), collapse = ", "))
  }
  tree <- read.tree(tree_path)
  asr_states <- readRDS(asr_states_rds_path)
  asr_states <- lapply(asr_states, function(v) {
    nv <- suppressWarnings(as.numeric(v)); names(nv) <- names(v); nv
  })

  n <- nrow(events)
  subtree_size <- integer(n)
  is_singleton <- logical(n)
  is_reacq     <- logical(n)
  reverts      <- rep(NA, n)
  cache <- new.env()   # events at the same node share one subtree walk

  for (i in seq_len(n)) {
    ev  <- events[i, ]
    key <- as.character(ev$descendant_node)
    if (is.null(cache[[key]])) {
      cache[[key]] <- get_subtree_within_distance(tree, ev$descendant_node,
                                                  max_branch_dist)
    }
    nodes <- cache[[key]]
    state_vec <- asr_states[[ev$character]]

    subtree_size[i] <- length(nodes)
    is_singleton[i] <- length(nodes) == 1

    path <- root_to_node_path(tree, ev$ancestral_node)
    is_reacq[i] <- count_state_changes(state_vec, path)$upstream_1_to_0_count > 0

    if (!is_singleton[i]) reverts[i] <- any_reversion_in_subtree(state_vec, nodes)
    if (i %% 500 == 0) message(i, " / ", n, " events")
  }

  out <- events[, need]
  out$subtree_size     <- subtree_size
  out$is_singleton     <- is_singleton
  out$is_reacquisition <- is_reacq
  out$subtree_reverts  <- reverts
  out <- out[, V2_COLUMNS]

  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(out, out_csv_path, row.names = FALSE)
  message("Singletons: ", sum(is_singleton), " (", round(100 * mean(is_singleton), 1),
          "%) | cluster members: ", sum(!is_singleton))
  invisible(out)
}

# Compare a rebuilt table with the shipped v2 table. All four should be 0.
compare_with_existing <- function(new_csv, old_csv) {
  new_tbl <- read.csv(new_csv, stringsAsFactors = FALSE)
  old_tbl <- read.csv(old_csv, stringsAsFactors = FALSE)
  stopifnot(nrow(new_tbl) == nrow(old_tbl))
  key_new <- paste(new_tbl$character, new_tbl$descendant_node)
  key_old <- paste(old_tbl$character, old_tbl$descendant_node)
  old_tbl <- old_tbl[match(key_new, key_old), ]
  for (col in c("subtree_size", "is_singleton", "is_reacquisition", "subtree_reverts")) {
    x <- new_tbl[[col]]; y <- old_tbl[[col]]
    mism <- sum((is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & x != y))
    cat(sprintf("%-18s mismatches vs shipped v2: %d\n", col, mism))
  }
}

# Usage (from the folder that contains phase4_out/):
#   source("scripts/tree_flag_acquisition_events.R")
#   run_final_exclusions()
#   compare_with_existing("phase4_out/v2_rebuilt.csv",
#                         "phase4_out/candidate_acquisition_events_flagged_v2.csv")