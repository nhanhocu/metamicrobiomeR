# Hierarchical FDR across taxonomic ranks

#' Parent of each feature in the taxonomy
#'
#' For a study whose features sit at several ranks (a `rank` column), finds for each
#' feature the feature at the nearest higher rank with the same lineage.
#'
#' @param study An `mm_study` with rank columns and a `rank` column.
#' @return A named character vector: the parent feature ID of each feature, `NA` for
#'   features without a parent in the study.
#' @family taxonomy
#' @export
feature_parents <- function(study) {
  .check_class(study, "mm_study", "study", "as_study()")
  f <- mm_features(study)
  if (!"rank" %in% colnames(f)) stop("The feature data has no 'rank' column.", call. = FALSE)
  cols <- intersect(.mm_ranks, colnames(f))
  path <- function(i, upto) {
    v <- unlist(f[i, cols[seq_len(upto)], drop = FALSE], use.names = FALSE)
    if (anyNA(v)) NA_character_ else paste(v, collapse = ";")
  }
  ri <- match(f$rank, cols)
  key <- vapply(seq_len(nrow(f)), function(i) if (is.na(ri[i])) NA_character_ else path(i, ri[i]), character(1))
  lookup <- stats::setNames(rownames(f), paste(f$rank, key, sep = "\r"))
  res <- vapply(seq_len(nrow(f)), function(i) {
    if (is.na(ri[i]) || ri[i] < 2) return(NA_character_)
    for (j in (ri[i] - 1):1) {
      kj <- path(i, j)
      hit <- if (is.na(kj)) NULL else lookup[paste(cols[j], kj, sep = "\r")]
      if (!is.null(hit) && !is.na(hit)) return(unname(hit))
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
  stats::setNames(res, rownames(f))
}

#' Hierarchical false discovery rate across ranks
#'
#' Tests the taxonomy top-down: the highest-level features are tested first with
#' Benjamini-Hochberg at level `q`; only the children of rejected features are tested
#' next, at a level scaled by the cumulative proportion of features rejected at the
#' higher levels (the selective adjustment of Benjamini and Bogomolov, applied level by
#' level as in tree-structured procedures such as TreeBH; the proportion is taken over
#' all features of a level, tested or not). Features whose parent was not rejected are not
#' tested. This replaces adjusting all ranks together, which ignores that a genus belongs
#' to a phylum. The error-control guarantee of the published procedures depends on
#' assumptions about the dependence between tests; check it by simulation for your design.
#'
#' @param pooled An `mm_pooled` object (see [pool_effects()]) that covers several ranks.
#' @param parent Named character vector mapping each feature ID to its parent ID (see
#'   [feature_parents()]). A feature whose ancestors are not in `pooled` is treated as a root.
#' @param q Target false discovery rate.
#' @param p_col Name of the p-value column.
#' @return `pooled` with the extra columns `level` (depth in the tree), `p_hier` (the
#'   hierarchically adjusted p-value, `NA` if the feature was not tested because its parent was
#'   not rejected) and `hier_rejected` (`p_hier <= q`).
#' @family pooling
#' @export
adjust_hierarchical <- function(pooled, parent, q = 0.05, p_col = "p_value") {
  .check_class(pooled, "mm_pooled", "pooled", "pool_effects()")
  if (!is.numeric(q) || q <= 0 || q >= 1) stop("'q' must be between 0 and 1.", call. = FALSE)
  if (is.null(names(parent))) stop("'parent' must be a named vector mapping each feature ID to its parent ID; feature_parents() makes one.", call. = FALSE)
  out <- as.data.frame(pooled)
  out$level <- NA_integer_; out$p_hier <- NA_real_; out$hier_rejected <- FALSE
  fam <- paste(out$part, out$term, sep = "\r")
  for (g in unique(fam)) {
    ii <- which(fam == g); feats <- out$feature[ii]
    anc <- function(f) {                                   # nearest ancestor present in this family
      p <- parent[f]
      while (!is.na(p) && !p %in% feats) p <- parent[p]
      if (is.na(p)) NA_character_ else unname(p)
    }
    par <- vapply(feats, anc, character(1), USE.NAMES = FALSE)
    depth <- integer(length(feats)); names(depth) <- feats
    for (k in seq_along(feats)) { d <- 1L; p <- par[k]; while (!is.na(p)) { d <- d + 1L; p <- par[match(p, feats)] }; depth[k] <- d }
    out$level[ii] <- depth
    pv <- out[[p_col]][ii]; names(pv) <- feats
    cum <- 1; rejected_prev <- character(0)
    for (lev in seq_len(max(depth))) {
      cand <- feats[depth == lev]
      n_level <- length(cand)                                # all features of the level, tested or not
      if (lev > 1) cand <- cand[par[match(cand, feats)] %in% rejected_prev]
      if (!length(cand)) break
      # BH within each family of siblings (children of one parent; all roots form one family)
      sib <- if (lev == 1) rep("root", length(cand)) else par[match(cand, feats)]
      rej <- character(0)
      for (s in unique(sib)) {
        cs <- cand[sib == s]; padj <- stats::p.adjust(pv[cs], "fdr") / cum
        out$p_hier[ii[match(cs, feats)]] <- pmin(padj, 1)
        rej <- c(rej, cs[padj <= q])
      }
      cum <- cum * length(rej) / n_level; rejected_prev <- rej
      if (!length(rej)) break
    }
  }
  out$hier_rejected <- !is.na(out$p_hier) & out$p_hier <= q
  .inform(sum(!is.na(out$p_hier)), " of ", nrow(out), " features were tested; ", sum(is.na(out$p_hier)),
          " were not because their parent was not rejected. ", sum(out$hier_rejected), " rejected at q = ", q, ".")
  class(out) <- c("mm_pooled", "data.frame")
  attr(out, "settings") <- c(attr(pooled, "settings"), list(hierarchical_q = q))
  out
}
