# Tests for R/12-taxonomy.R
test_that("parse_taxonomy understands the common conventions", {
  gg <- parse_taxonomy("k__Bacteria; p__Firmicutes; c__Clostridia; o__; f__; g__; s__")
  expect_equal(gg$kingdom, "Bacteria"); expect_equal(gg$class, "Clostridia")
  expect_true(all(is.na(gg[, c("order", "family", "genus", "species")])))
  expect_equal(gg$deepest, "class")
  silva <- parse_taxonomy("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella;s__uncultured_bacterium")
  expect_equal(silva$genus, "Escherichia-Shigella"); expect_true(is.na(silva$species)); expect_equal(silva$deepest, "genus")
  old <- parse_taxonomy("D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales")
  expect_equal(old$order, "Lactobacillales"); expect_equal(old$deepest, "order")
  dada <- parse_taxonomy("Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;NA")
  expect_equal(dada$genus, "Bacteroides"); expect_true(is.na(dada$species))
  mp <- parse_taxonomy("k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Blautia|s__Blautia_obeum|t__SGB1")
  expect_equal(mp$species, "Blautia_obeum"); expect_equal(mp$deepest, "species")
  un <- parse_taxonomy(c("Unassigned", NA, "", "k__Bacteria;p__Firmicutes_unclassified"))
  expect_true(all(is.na(un$deepest[1:3]))); expect_true(is.na(un$phylum[4])); expect_equal(un$deepest[4], "kingdom")
  expect_s3_class(parse_taxonomy("k__A"), "data.frame")
})

test_that("harmonise_names removes database decorations and applies synonyms", {
  expect_equal(harmonise_names(c("[Ruminococcus]", "Firmicutes_A", "Escherichia Shigella", "Escherichia-Shigella")),
               c("Ruminococcus", "Firmicutes", "Escherichia_Shigella", "Escherichia_Shigella"))
  expect_equal(harmonise_names(c("Bacillota", "bacteroidota", "Pseudomonadota")), c("Firmicutes", "Bacteroidetes", "Proteobacteria"))
  expect_equal(harmonise_names("Foo", synonyms = c(Foo = "Bar")), "Bar")
  expect_true(is.na(harmonise_names(NA)))
})

toy_asv <- function() {
  set.seed(4)
  counts <- matrix(rpois(7 * 6, 30), 7, 6, dimnames = list(paste0("ASV", 1:7), paste0("S", 1:6)))
  tax <- data.frame(kingdom = "Bacteria",
                    phylum = c("Firmicutes", "Firmicutes", "Firmicutes", "Bacteroidota", "Bacteroidota", "Firmicutes_A", NA),
                    genus  = c("Blautia", "Blautia", "[Ruminococcus]", "Bacteroides", NA, "Ruminococcus", NA),
                    row.names = rownames(counts), stringsAsFactors = FALSE)
  as_study(counts = counts, taxonomy = tax, samples = data.frame(g = 1:6, row.names = colnames(counts)),
           study = list(id = "toy", feature_type = "asv"))
}

test_that("agglomerate sums features per taxon and keeps unassigned reads", {
  s <- toy_asv(); a <- agglomerate(s, "genus")
  expect_s3_class(a, "mm_study")
  expect_equal(unname(colSums(mm_counts(a))), unname(colSums(mm_counts(s))))            # totals preserved
  expect_equal(unname(mm_counts(a)["Bacteria;Firmicutes;Blautia", ]), unname(colSums(mm_counts(s)[1:2, ])))
  expect_true("unassigned_genus" %in% rownames(mm_counts(a)))
  expect_equal(unname(mm_counts(a)["unassigned_genus", ]), unname(colSums(mm_counts(s)[c(5, 7), ])))
  f <- mm_features(a)
  expect_true(all(f$rank == "genus"))
  expect_true(is.na(f["unassigned_genus", "genus"]))
  expect_equal(unname(mm_depth(a)), unname(mm_depth(s)))
  d <- agglomerate(s, "genus", unassigned = "drop")
  expect_false("unassigned_genus" %in% rownames(mm_counts(d)))
  expect_lt(sum(mm_counts(d)), sum(mm_counts(s)))
})

test_that("key = 'name' merges equal names across lineages and harmonises them", {
  s <- toy_asv(); a <- agglomerate(s, "genus", key = "name")
  expect_true("Ruminococcus" %in% rownames(mm_counts(a)))
  expect_equal(unname(mm_counts(a)["Ruminococcus", ]), unname(colSums(mm_counts(s)[c(3, 6), ])))   # [Ruminococcus] + Ruminococcus
  p <- agglomerate(s, "phylum", key = "name")
  expect_equal(unname(mm_counts(p)["Firmicutes", ]), unname(colSums(mm_counts(s)[c(1, 2, 3, 6), ])))   # Firmicutes_A merged
})

test_that("agglomeration works on proportions and fit_study then leaves out the unassigned feature", {
  s <- toy_asv()
  rel <- as_study(relabund = mm_relabund(s), taxonomy = mm_features(s), samples = mm_samples(s), study = list(id = "r", feature_type = "asv"))
  a <- agglomerate(rel, "phylum")
  expect_equal(unname(colSums(mm_relabund(a))), rep(1, 6), tolerance = 1e-10)
  expect_null(mm_counts(a))
  smp <- cbind(mm_samples(s), grp = factor(rep(c("a", "b"), 3)))
  cnt <- as_study(counts = mm_counts(s), taxonomy = mm_features(s), samples = smp, study = list(id = "c", feature_type = "asv"))
  e <- fit_study(agglomerate(cnt, "genus"), ~ grp, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0)
  expect_false("unassigned_genus" %in% e$feature)
})

test_that("agglomerate validates its input", {
  s <- toy_asv()
  expect_error(agglomerate(list(), "genus"), "mm_study")
  expect_error(agglomerate(s, "family"), "no 'family' column")
  expect_error(agglomerate(s, "nonsense"))
  expect_error(agglomerate(agglomerate(s, "genus"), "phylum"), "already agglomerated")
})

# ---- the vectorised parser must equal the original loop implementation
parse_taxonomy_reference <- function(x, sep = NULL) {
  x <- as.character(x)
  if (is.null(sep)) sep <- if (any(grepl("|", x, fixed = TRUE))) "|" else ";"
  ranks <- metamicrobiomeR:::.mm_ranks; pr <- metamicrobiomeR:::.mm_prefix_rank
  out <- matrix(NA_character_, length(x), length(ranks), dimnames = list(NULL, ranks))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || !nzchar(trimws(x[i]))) next
    parts <- trimws(strsplit(x[i], sep, fixed = TRUE)[[1]])
    silva_old <- grepl("^D_[0-9]+__", parts); prefixed <- grepl("^[A-Za-z]__", parts)
    if (all(silva_old)) { pos <- as.integer(sub("^D_([0-9]+)__.*", "\\1", parts)) + 1L; val <- sub("^D_[0-9]+__", "", parts)
    } else if (all(prefixed)) { rk <- match(tolower(substr(parts, 1, 1)), names(pr)); keep <- !is.na(rk); pos <- match(pr[rk[keep]], ranks); val <- substring(parts[keep], 4)
    } else { pos <- seq_along(parts); val <- parts }
    ok <- pos >= 1 & pos <= length(ranks)
    out[i, pos[ok]] <- metamicrobiomeR:::.clean_taxon(val[ok])
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$deepest <- vapply(seq_len(nrow(out)), function(i) { k <- which(!is.na(unlist(out[i, ranks], use.names = FALSE))); if (length(k)) ranks[max(k)] else NA_character_ }, character(1))
  out
}

test_that("the vectorised parse_taxonomy equals the original implementation on all conventions", {
  set.seed(42)
  nm <- function(n) sample(c("Firmicutes", "Bacteroidota", "[Ruminococcus]", "Escherichia-Shigella", "uncultured", "Unassigned", "unclassified", "other", "", "Clostridia_unclassified",
                             "Lachnospiraceae", "Blautia", "Faecalibacterium_prausnitzii", "uncultured_bacterium", "NA"), n, TRUE)
  gen <- function(style, n) vapply(seq_len(n), function(i) {
    len <- sample(1:7, 1); v <- nm(len)
    switch(style,
      greengenes = paste(paste0(c("k", "p", "c", "o", "f", "g", "s")[seq_len(len)], "__", v), collapse = "; "),
      silva_q2 = paste(paste0(c("d", "p", "c", "o", "f", "g", "s")[seq_len(len)], "__", v), collapse = ";"),
      silva_old = paste(paste0("D_", seq_len(len) - 1, "__", v), collapse = ";"),
      positional = paste(v, collapse = ";"),
      mixed = paste(c(paste0("k__", v[1]), v[-1]), collapse = ";"),
      metaphlan = paste(c(paste0(c("k", "p", "c", "o", "f", "g", "s")[seq_len(len)], "__", v), if (runif(1) < 0.3) "t__SGB1"), collapse = "|"))
  }, character(1))
  for (style in c("greengenes", "silva_q2", "silva_old", "positional", "mixed")) {
    x <- c(gen(style, 400), NA, "", "   ", "Unassigned")
    expect_identical(parse_taxonomy(x), parse_taxonomy_reference(x), info = style)
  }
  xm <- c(gen("metaphlan", 400), "k__Bacteria")
  expect_identical(parse_taxonomy(xm), parse_taxonomy_reference(xm), info = "metaphlan")
  xs <- gen("greengenes", 50); xs <- gsub(";", "|", xs)
  expect_identical(parse_taxonomy(xs, sep = "|"), parse_taxonomy_reference(xs, sep = "|"))
  expect_identical(parse_taxonomy(character(0)), parse_taxonomy_reference(character(0)))
})
