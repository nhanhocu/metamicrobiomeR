# Taxonomy: parse the naming conventions of different databases, harmonise names, find parents

.mm_prefix_rank <- c(d = "kingdom", k = "kingdom", p = "phylum", c = "class", o = "order",
                     f = "family", g = "genus", s = "species")

.mm_placeholder <- "^(|na|other|unassigned|unclassified|unidentified|uncultured|unknown|ambiguous_taxa|metagenome)$|(_|\\s)(unclassified|unassigned)$|^(uncultured|unidentified)(_|\\s)"

.clean_taxon <- function(v) {
  v <- trimws(v)
  v[is.na(v) | grepl(.mm_placeholder, tolower(v))] <- NA_character_
  v
}

#' Parse lineage strings into rank columns
#'
#' Understands the common conventions: prefixed ranks (`k__`/`d__`, `p__`, `c__`, `o__`,
#' `f__`, `g__`, `s__` as in Greengenes, SILVA via QIIME 2, GTDB and MetaPhlAn), the older
#' SILVA style `D_0__`, `D_1__`, ..., and plain lineages whose ranks are given by position
#' (as from DADA2's `assignTaxonomy`). Placeholders for unassigned levels (empty names,
#' `Unassigned`, `uncultured`, `unclassified`, `other`, ...) become `NA`.
#'
#' @param x Character vector of lineage strings.
#' @param sep Separator between ranks; by default `"|"` if any string contains it, else `";"`.
#' @return A data.frame with the columns kingdom, phylum, class, order, family, genus and
#'   species (character, `NA` where unassigned), and `deepest`, the deepest rank that is assigned.
#' @examples
#' silva <- paste("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria",
#'                "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella", sep = ";")
#' parse_taxonomy(c("k__Bacteria; p__Firmicutes; c__Clostridia; o__; f__; g__", silva, "Unassigned"))
#' @family taxonomy
#' @export
parse_taxonomy <- function(x, sep = NULL) {
  x <- as.character(x)
  if (is.null(sep)) sep <- if (any(grepl("|", x, fixed = TRUE))) "|" else ";"
  n <- length(x); nr <- length(.mm_ranks)
  out <- matrix(NA_character_, n, nr, dimnames = list(NULL, .mm_ranks))
  ok_str <- !is.na(x) & nzchar(trimws(x))
  parts <- strsplit(x[ok_str], sep, fixed = TRUE)
  lens <- lengths(parts)
  if (length(parts) && any(lens > 0)) {
    tok <- trimws(unlist(parts)); id <- rep(which(ok_str), lens); at <- sequence(lens)
    is_old <- grepl("^D_[0-9]+__", tok); is_pref <- grepl("^[A-Za-z]__", tok)
    str_old <- tapply(is_old, id, all)[as.character(id)]; str_pref <- tapply(is_pref, id, all)[as.character(id)]
    pos <- at; val <- tok
    o <- str_old; pos[o] <- as.integer(sub("^D_([0-9]+)__.*", "\\1", tok[o])) + 1L; val[o] <- sub("^D_[0-9]+__", "", tok[o])
    pf <- !str_old & str_pref
    rk <- match(tolower(substr(tok[pf], 1, 1)), names(.mm_prefix_rank))
    pos[pf] <- match(.mm_prefix_rank[rk], .mm_ranks); pos[pf][is.na(rk)] <- NA_integer_       # e.g. MetaPhlAn strain level t__ is ignored
    val[pf] <- substring(tok[pf], 4)
    keep <- !is.na(pos) & pos >= 1 & pos <= nr
    out[cbind(id[keep], pos[keep])] <- .clean_taxon(val[keep])
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  filled <- !is.na(as.matrix(out))
  out$deepest <- ifelse(rowSums(filled) > 0, .mm_ranks[pmax.int(1L, max.col(filled * rep(seq_len(nr), each = n), ties.method = "last"))], NA_character_)
  out$deepest[rowSums(filled) == 0] <- NA_character_
  out
}

#' Harmonise taxon names across databases
#'
#' Removes the decorations that differ between taxonomies so that the same genus or
#' phylum can be matched across studies: square brackets (`[Ruminococcus]`), GTDB
#' suffixes (`Firmicutes_A`), underscores versus spaces, and capitalisation; then
#' applies a synonym table. The default table maps the phylum names that were
#' renamed in recent taxonomies to the older ones (`Bacillota` to `Firmicutes`, ...).
#' Harmonisation across GTDB and SILVA or Greengenes is lossy: GTDB splits some taxa
#' (for example Firmicutes into several phyla) and renames many genera, so check the
#' matched names for the ranks you analyse.
#'
#' @param x Character vector of taxon names.
#' @param synonyms Named character vector (names are the alternative names, values the
#'   names to use) added to the default table, which it overrides. `NULL` uses only the default.
#' @return Character vector of the same length (`NA` stays `NA`).
#' @examples
#' harmonise_names(c("[Ruminococcus]", "Firmicutes_A", "Bacillota", "Escherichia Shigella"))
#' @family taxonomy
#' @export
harmonise_names <- function(x, synonyms = NULL) {
  syn <- c(Bacillota = "Firmicutes", Actinomycetota = "Actinobacteria", Pseudomonadota = "Proteobacteria",
           Bacteroidota = "Bacteroidetes", Fusobacteriota = "Fusobacteria", Verrucomicrobiota = "Verrucomicrobia",
           Mycoplasmatota = "Tenericutes", Campylobacterota = "Epsilonbacteraeota")
  if (!is.null(synonyms)) syn[names(synonyms)] <- synonyms
  y <- as.character(x)
  y <- gsub("^\\[|\\]$|\\[|\\]", "", y)
  y <- sub("_[A-Z]{1,2}$", "", y)                          # GTDB suffix (Firmicutes_A)
  y <- gsub("[ ]+", "_", trimws(y))
  y <- gsub("-", "_", y)
  hit <- match(tolower(y), tolower(names(syn)))
  y[!is.na(hit)] <- unname(syn)[hit[!is.na(hit)]]
  y[is.na(x)] <- NA_character_
  y
}

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
