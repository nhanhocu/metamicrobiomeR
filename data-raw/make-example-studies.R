# Builds data/example_studies.rda: four SIMULATED gut-microbiome studies on different platforms.
# Nothing here is real data. The genus names are real, the abundances, effects and covariates are invented.
suppressMessages(library(metamicrobiomeR))
set.seed(20261)
tax <- read.table(text = "
genus phylum
Bacteroides Bacteroidetes
Prevotella Bacteroidetes
Alistipes Bacteroidetes
Parabacteroides Bacteroidetes
Odoribacter Bacteroidetes
Barnesiella Bacteroidetes
Faecalibacterium Firmicutes
Roseburia Firmicutes
Blautia Firmicutes
Ruminococcus Firmicutes
Coprococcus Firmicutes
Dorea Firmicutes
Lachnospira Firmicutes
Eubacterium Firmicutes
Anaerostipes Firmicutes
Butyricicoccus Firmicutes
Subdoligranulum Firmicutes
Oscillibacter Firmicutes
Gemmiger Firmicutes
Flavonifractor Firmicutes
Lachnoclostridium Firmicutes
Clostridium Firmicutes
Romboutsia Firmicutes
Intestinibacter Firmicutes
Turicibacter Firmicutes
Holdemanella Firmicutes
Catenibacterium Firmicutes
Christensenella Firmicutes
Veillonella Firmicutes
Dialister Firmicutes
Megasphaera Firmicutes
Phascolarctobacterium Firmicutes
Acidaminococcus Firmicutes
Enterococcus Firmicutes
Streptococcus Firmicutes
Lactobacillus Firmicutes
Lactococcus Firmicutes
Gemella Firmicutes
Granulicatella Firmicutes
Staphylococcus Firmicutes
Bifidobacterium Actinobacteria
Collinsella Actinobacteria
Eggerthella Actinobacteria
Rothia Actinobacteria
Actinomyces Actinobacteria
Corynebacterium Actinobacteria
Escherichia Proteobacteria
Klebsiella Proteobacteria
Citrobacter Proteobacteria
Enterobacter Proteobacteria
Proteus Proteobacteria
Pseudomonas Proteobacteria
Haemophilus Proteobacteria
Sutterella Proteobacteria
Desulfovibrio Proteobacteria
Bilophila Proteobacteria
Akkermansia Verrucomicrobia
Fusobacterium Fusobacteria
", header = TRUE, stringsAsFactors = FALSE)
G <- tax$genus; ng <- length(G)
# invented effect of 'case' versus 'control' on the log-odds of each genus
effect <- setNames(rep(0, ng), G)
effect[c("Escherichia", "Fusobacterium", "Enterococcus")] <- 1.5        # higher in cases
effect[c("Faecalibacterium", "Roseburia", "Butyricicoccus")] <- -1.3   # lower in cases
base <- setNames(rnorm(ng, -5, 1.4), G); base[c("Bacteroides", "Faecalibacterium", "Blautia", "Prevotella")] <- c(-2.2, -2.6, -2.8, -3.0)

make_study <- function(id, genera, depth, platform, region, db, asv, rename, n = 50, bias_sd = 0.3, effect_scale = 1, person_sd = 1.0) {
  bias <- setNames(rnorm(length(genera), 0, bias_sd), genera)             # platform-specific distortion of each genus
  group <- factor(rep(c("control", "case"), each = n / 2), levels = c("control", "case"))
  age <- round(runif(n, 25, 75)); sex <- factor(sample(c("F", "M"), n, TRUE))
  # genus-level expected proportions per sample; each genus then split over a few ASVs for the 16S studies
  nasv <- if (asv) sample(1:4, length(genera), TRUE) else rep(1L, length(genera))
  feat_genus <- rep(genera, nasv); feat_id <- if (asv) sprintf("%s_ASV%03d", id, seq_along(feat_genus)) else paste0(id, "_", feat_genus)
  split <- unlist(lapply(nasv, function(k) { w <- rgamma(k, 2); w / sum(w) }))
  counts <- matrix(0L, length(feat_genus), n, dimnames = list(feat_id, sprintf("%s_S%02d", id, seq_len(n))))
  for (j in seq_len(n)) {
    lp <- base[genera] + bias + rnorm(length(genera), 0, person_sd) +          # person-to-person variation in every genus
      (group[j] == "case") * effect_scale * effect[genera] + 0.004 * (age[j] - 50) * (genera == "Bacteroides")
    p <- exp(lp) / sum(exp(lp)); pf <- p[feat_genus] * split
    counts[, j] <- rmultinom(1, round(depth * exp(rnorm(1, 0, 0.35))), pf)
  }
  as_study(counts = counts, taxonomy = data.frame(phylum = rename$phylum(feat_genus), genus = rename$genus(feat_genus), row.names = feat_id, stringsAsFactors = FALSE),
           samples = data.frame(group = group, age = age, sex = sex, row.names = colnames(counts)),
           study = list(id = id, feature_type = if (asv) "asv" else "taxon", platform = platform, region = region, database = db))
}
ph <- function(g) tax$phylum[match(g, tax$genus)]
keep <- function(a, b) G[a:b]
example_studies <- list(
  CohortA = make_study("CohortA", keep(1, 52), 3e4, "Illumina MiSeq", "16S V4", "SILVA 138", TRUE,
                       list(phylum = ph, genus = function(g) ifelse(g == "Eubacterium", "[Eubacterium]", ifelse(g == "Ruminococcus", "[Ruminococcus]", g)))),
  CohortB = make_study("CohortB", G[c(1:12, 14:57)], 2e4, "Ion Torrent", "16S V3-V4", "SILVA 138", TRUE, list(phylum = ph, genus = identity), bias_sd = 0.4),
  CohortC = make_study("CohortC", G[c(1:10, 13:58)], 5e5, "NovaSeq shotgun", NA_character_, "GTDB r207", FALSE,
                       list(phylum = function(g) ifelse(ph(g) == "Firmicutes", "Firmicutes_A", ph(g)),
                            genus = function(g) ifelse(g %in% c("Blautia", "Ruminococcus", "Clostridium"), paste0(g, "_A"), g)), bias_sd = 0.3))
example_studies$CohortD <- make_study("CohortD", G[2:55], 4e4, "Illumina MiSeq", "16S V1-V2", "SILVA 138", TRUE, list(phylum = ph, genus = identity), bias_sd = 0.3)
for (s in example_studies) check_study(s)
dir.create("data", showWarnings = FALSE)
save(example_studies, file = "data/example_studies.rda", compress = "xz")
cat("saved", format(file.size("data/example_studies.rda") / 1024, digits = 3), "KB\n")
