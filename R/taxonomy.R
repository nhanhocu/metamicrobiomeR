# Taxonomy: parsing the naming conventions of different databases, harmonising names,
# and agglomerating features to a taxonomic rank

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
#' parse_taxonomy(c("k__Bacteria; p__Firmicutes; c__Clostridia; o__; f__; g__",
#'                  "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella",
#'                  "Unassigned"))
#' @export
parse_taxonomy <- function(x, sep = NULL) {
  x <- as.character(x)
  if (is.null(sep)) sep <- if (any(grepl("|", x, fixed = TRUE))) "|" else ";"
  out <- matrix(NA_character_, length(x), length(.mm_ranks), dimnames = list(NULL, .mm_ranks))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || !nzchar(trimws(x[i]))) next
    parts <- trimws(strsplit(x[i], sep, fixed = TRUE)[[1]])
    silva_old <- grepl("^D_[0-9]+__", parts)
    prefixed <- grepl("^[A-Za-z]__", parts)
    if (all(silva_old)) {
      pos <- as.integer(sub("^D_([0-9]+)__.*", "\\1", parts)) + 1L
      val <- sub("^D_[0-9]+__", "", parts)
    } else if (all(prefixed)) {
      rk <- match(tolower(substr(parts, 1, 1)), names(.mm_prefix_rank))
      keep <- !is.na(rk)                                    # e.g. MetaPhlAn strain level t__ is ignored
      pos <- match(.mm_prefix_rank[rk[keep]], .mm_ranks); val <- substring(parts[keep], 4)
    } else {
      pos <- seq_along(parts); val <- parts
    }
    ok <- pos >= 1 & pos <= length(.mm_ranks)
    out[i, pos[ok]] <- .clean_taxon(val[ok])
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$deepest <- vapply(seq_len(nrow(out)), function(i) {
    k <- which(!is.na(unlist(out[i, .mm_ranks], use.names = FALSE)))
    if (length(k)) .mm_ranks[max(k)] else NA_character_
  }, character(1))
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

#' Agglomerate the features of a study to a taxonomic rank
#'
#' Sums the abundances of features that share a taxon at the requested rank (for
#' example ASVs into genera), using the rank columns of the feature data. Features
#' without an assignment at that rank are collected into one feature
#' (`unassigned_<rank>`) so that sample totals, and therefore proportions and
#' sequencing depth, stay unchanged; [fit_study()] leaves such features out of the
#' models by default.
#'
#' @param study An `mm_study` with taxonomy in rank columns whose features are not
#'   yet agglomerated (no `rank` column, or `rank = "leaf"`).
#' @param rank One of kingdom, phylum, class, order, family, genus, species.
#' @param key How taxa are identified: `"lineage"` (default) uses the lineage down to the
#'   rank, so equal names in different lineages stay separate; `"name"` uses the (harmonised)
#'   name at the rank only, which is what matching across studies and databases needs.
#' @param unassigned `"keep"` (default) or `"drop"` the unassigned feature. Dropping
#'   changes the sample totals and so the proportions.
#' @param synonyms Passed to [harmonise_names()] when `key = "name"`.
#' @return An `mm_study` with one feature per taxon, a `rank` column, and the rank
#'   columns down to the requested rank.
#' @export
agglomerate <- function(study, rank, key = c("lineage", "name"), unassigned = c("keep", "drop"), synonyms = NULL) {
  if (!inherits(study, "mm_study")) stop("'study' must be an mm_study object.", call. = FALSE)
  key <- match.arg(key); unassigned <- match.arg(unassigned)
  rank <- match.arg(rank, .mm_ranks)
  feats <- mm_features(study)
  if (!rank %in% colnames(feats)) stop("The feature data has no '", rank, "' column.", call. = FALSE)
  if ("rank" %in% colnames(feats) && any(!is.na(feats$rank) & feats$rank != "leaf")) {
    stop("The study is already agglomerated (its features have ranks); choose a rank in fit_study() instead.", call. = FALSE)
  }
  cols <- intersect(.mm_ranks[seq_len(match(rank, .mm_ranks))], colnames(feats))
  tax <- feats[, cols, drop = FALSE]
  nm <- tax[[rank]]
  id <- if (key == "name") harmonise_names(nm, synonyms) else {
    ifelse(is.na(nm), NA_character_, apply(tax, 1, function(r) paste(r[!is.na(r)], collapse = ";")))
  }
  un <- is.na(id)
  id[un] <- paste0("unassigned_", rank)
  if (unassigned == "drop") {
    keep <- !un
    if (!any(keep)) stop("No feature is assigned at rank '", rank, "'.", call. = FALSE)
  } else keep <- rep(TRUE, length(id))
  agg <- function(m) if (is.null(m)) NULL else rowsum(m[keep, , drop = FALSE], group = id[keep], reorder = FALSE)
  counts <- agg(mm_counts(study)); rel <- agg(study$assays$relabund)
  first <- !duplicated(id[keep])
  newf <- tax[keep, , drop = FALSE][first, , drop = FALSE]
  rownames(newf) <- id[keep][first]
  if (key == "name") newf[, setdiff(cols, rank)] <- NA_character_
  newf[[rank]] <- if (key == "name") rownames(newf) else newf[[rank]]
  newf[grepl("^unassigned_", rownames(newf)), cols] <- NA_character_
  newf$rank <- rank
  info <- study$study; info$feature_type <- "taxon"
  as_study(counts = counts, relabund = rel, taxonomy = newf, samples = mm_samples(study), study = info)
}
