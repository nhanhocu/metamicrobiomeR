# Simulated studies and the example data

.softmax <- function(x) { e <- exp(x - max(x)); e / sum(e) }

.sim_one <- function(base, beta, n_samples, depth, depth_sd, overdispersion, zero_inflation, n_subjects, subject_sd, between, id) {
  nf <- length(base)
  person <- if (is.null(n_subjects)) NULL else rep(seq_len(n_subjects), length.out = n_samples)
  if (!is.null(person)) person <- sort(person)
  grp <- if (!is.null(person) && between) factor(ifelse(person <= n_subjects / 2, "a", "b")) else factor(rep(c("a", "b"), length.out = n_samples))
  off <- if (is.null(person)) NULL else matrix(stats::rnorm(n_subjects * nf, 0, subject_sd), n_subjects, nf)
  N <- pmax(round(exp(stats::rnorm(n_samples, log(depth), depth_sd))), 100)
  cnt <- matrix(0L, nf, n_samples, dimnames = list(paste0("F", seq_len(nf)), paste0("S", seq_len(n_samples))))
  for (i in seq_len(n_samples)) {
    lp <- base + beta * (grp[i] == "b") + if (is.null(off)) 0 else off[person[i], ]
    p <- .softmax(lp)
    g <- stats::rgamma(nf, shape = pmax(overdispersion * p, 1e-3)); p <- g / sum(g)
    y <- as.integer(stats::rmultinom(1, N[i], p))
    if (zero_inflation > 0) y[stats::rbinom(nf, 1, zero_inflation) == 1] <- 0L
    cnt[, i] <- y
  }
  smp <- data.frame(group = grp, row.names = colnames(cnt))
  if (!is.null(person)) smp$person <- person
  as_study(counts = cnt, samples = smp, study = list(id = id, feature_type = "asv", subject = if (!is.null(person)) "person"))
}

#' Simulate a study with known group differences
#'
#' Draws counts from a Dirichlet-multinomial model: baseline proportions on a log scale, a
#' group effect on the log proportion of some features, sample-specific overdispersion,
#' variable sequencing depth, optional structural zeros and optional repeated measures with a
#' subject effect. Because the data are compositional, a change in some features shifts the
#' proportions of all the others; with `balanced = FALSE` (all differential features move in
#' the same direction) this compositional effect is strong.
#'
#' @param n_samples,n_features Number of samples and features.
#' @param n_diff Number of truly differential features.
#' @param effect Group effect on the log proportion of the differential features.
#' @param balanced If `TRUE`, differential features alternate between increasing and decreasing.
#' @param depth,depth_sd Median sequencing depth and the standard deviation of its logarithm.
#' @param overdispersion Dirichlet concentration (smaller means more overdispersion).
#' @param zero_inflation Probability that a count is set to zero.
#' @param n_subjects Number of subjects for repeated measures (`NULL` for independent samples).
#' @param subject_sd Standard deviation of the subject effect on the log proportions.
#' @param between If `TRUE` with subjects, the group is constant within a subject.
#' @param id Study ID.
#' @param seed Optional seed.
#' @return A list with the `mm_study` (`study`) and the `truth`: a data.frame with the feature,
#'   whether it is differential (`diff`) and its true log effect (`beta`).
#' @family simulation
#' @export
simulate_study <- function(n_samples = 60, n_features = 30, n_diff = 6, effect = 1, balanced = TRUE, depth = 20000,
                           depth_sd = 0.5, overdispersion = 20, zero_inflation = 0, n_subjects = NULL, subject_sd = 0.5,
                           between = FALSE, id = "sim", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (n_diff > n_features) stop("'n_diff' cannot exceed 'n_features'.", call. = FALSE)
  base <- sort(stats::rnorm(n_features, 0, 1.5), decreasing = TRUE)
  di <- sort(sample.int(n_features, n_diff))
  beta <- numeric(n_features)
  beta[di] <- effect * if (balanced) rep(c(1, -1), length.out = n_diff) else 1
  study <- .sim_one(base, beta, n_samples, depth, depth_sd, overdispersion, zero_inflation, n_subjects, subject_sd, between, id)
  list(study = study, truth = data.frame(feature = paste0("F", seq_len(n_features)), diff = beta != 0, beta = beta, stringsAsFactors = FALSE))
}

#' Simulate several studies with between-study heterogeneity
#'
#' Studies share the features, the baseline composition (up to a study-specific perturbation)
#' and the set of differential features, but each has its own effect size, drawn around
#' `effect` with standard deviation `tau`. Optionally the last study is a small, sparse study.
#'
#' @param n_studies Number of studies.
#' @param tau Between-study standard deviation of the effect.
#' @param small_study If `TRUE`, the last study has only a fifth of the samples.
#' @inheritParams simulate_study
#' @return A list with `studies` (a list of `mm_study` objects) and `truth` (as in [simulate_study()]).
#' @family simulation
#' @export
simulate_studies <- function(n_studies = 4, tau = 0.3, small_study = FALSE, n_samples = 60, n_features = 30, n_diff = 6,
                             effect = 1, balanced = TRUE, depth = 20000, depth_sd = 0.5, overdispersion = 20,
                             zero_inflation = 0, n_subjects = NULL, subject_sd = 0.5, between = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  base <- sort(stats::rnorm(n_features, 0, 1.5), decreasing = TRUE)
  di <- sort(sample.int(n_features, n_diff))
  sign_ <- if (balanced) rep(c(1, -1), length.out = n_diff) else rep(1, n_diff)
  studies <- vector("list", n_studies)
  for (k in seq_len(n_studies)) {
    beta <- numeric(n_features); beta[di] <- (effect + stats::rnorm(1, 0, tau)) * sign_
    n_k <- if (small_study && k == n_studies) max(round(n_samples / 5), 8) else n_samples
    studies[[k]] <- .sim_one(base + stats::rnorm(n_features, 0, 0.3), beta, n_k, depth, depth_sd, overdispersion, zero_inflation,
                             if (is.null(n_subjects)) NULL else min(n_subjects, n_k %/% 2), subject_sd, between, paste0("sim", k))
  }
  truth <- data.frame(feature = paste0("F", seq_len(n_features)), diff = FALSE, beta = 0, stringsAsFactors = FALSE)
  truth$diff[di] <- TRUE; truth$beta[di] <- effect * sign_
  list(studies = studies, truth = truth)
}

#' Simulated example studies on three platforms
#'
#' Four **simulated** gut-microbiome case-control studies for trying out the package. The genus
#' names are real; the abundances, effects and covariates are invented, so nothing here says
#' anything about real biology. They were made to look like a cross-platform meta-analysis:
#'
#' \describe{
#'   \item{CohortA}{16S V4, Illumina MiSeq, SILVA 138, ASV level (50 samples)}
#'   \item{CohortB}{16S V3-V4, Ion Torrent, SILVA 138, ASV level}
#'   \item{CohortC}{shotgun, NovaSeq, GTDB r207, genus level (names such as `Blautia_A`)}
#'   \item{CohortD}{16S V1-V2, Illumina MiSeq, SILVA 138, ASV level}
#' }
#'
#' Each study has counts, a taxonomy (phylum and genus), and the sample variables `group`
#' (control or case), `age` and `sex`. The studies cover different sets of genera, use different
#' naming conventions (`[Eubacterium]`, `Blautia_A`), very different sequencing depths and a
#' platform-specific distortion of each genus, as in real cross-platform data.
#'
#' Six genera were given a known difference between cases and controls: *Escherichia*,
#' *Enterococcus* and *Fusobacterium* are higher, *Faecalibacterium*, *Roseburia* and
#' *Butyricicoccus* lower. *Fusobacterium* occurs in one study only, so it cannot be pooled.
#' Every other genus has no true effect, although its relative abundance shifts a little because
#' relative abundances are compositional.
#'
#' @format A named list of four `mm_study` objects.
#' @source Simulated by `data-raw/make-example-studies.R`.
#' @examples
#' study_table(example_studies)[, c("study", "n_samples", "platform", "database")]
#' check_poolable(example_studies, rank = "genus")
#' res <- meta_analyse(example_studies, ~ group + age, term = "groupcase", rank = "genus",
#'                     engine = "lm_asin", fit_args = list(longitudinal = FALSE))
#' res
"example_studies"
