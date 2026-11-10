suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); fdir <- file.path(S, "crc", "fits"); odir <- file.path(S, "crc", "pooled")
st <- readRDS(file.path(S, "crc", "crc_studies_genus.rds"))
plat <- vapply(st, function(x) x$study$platform, character(1))
eng <- commandArgs(TRUE); if (!length(eng)) eng <- c("clr", "asin", "zibb", "clr_adj", "bezi", "zibb_adj")
pc <- c("Fusobacterium", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "Gemella")
load_eff <- function(e, studies = NULL) {
  fs <- list.files(fdir, paste0("^", e, "__"), full.names = TRUE)
  ids <- sub(paste0("^", e, "__(.*)\\.rds$"), "\\1", basename(fs))
  keep <- if (is.null(studies)) TRUE else ids %in% studies
  xs <- lapply(fs[keep], readRDS); xs <- xs[!vapply(xs, is.null, logical(1))]
  do.call(bind_effects, xs)
}
q <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"), message = function(m) invokeRestart("muffleMessage"))
summ <- list()
for (e in eng) {
  if (!length(list.files(fdir, paste0("^", e, "__")))) next
  sets <- list(all = NULL, shotgun = names(plat)[plat == "shotgun"], `16S` = setdiff(names(plat)[plat == "16S"], "ZellerG_2014_16S"))
  if (grepl("^zibb", e)) sets <- sets["16S"]           # the count model did not converge on MetaPhlAn pseudo-counts
  for (sn in names(sets)) {
    eff <- load_eff(e, sets[[sn]]); n_raw <- nrow(eff)
    effq <- q(filter_effects(eff)); effq <- structure(effq, class = class(eff))
    for (cfg in list(c("REML", "adhoc"), c("REML", "z"), c("DL", "z"))) {
      p <- q(pool_effects(effq, term = "groupCRC", method = cfg[1], test = cfg[2]))
      p <- as.data.frame(p); p$engine <- e; p$set <- sn; p$method <- cfg[1]; p$test <- cfg[2]
      tag <- paste(e, sn, cfg[1], cfg[2], sep = "_")
      saveRDS(p, file.path(odir, paste0(tag, ".rds")))
      pcs <- p[p$feature %in% pc, ]
      summ[[tag]] <- data.frame(engine = e, set = sn, method = cfg[1], test = cfg[2], studies = length(unique(eff$study)),
        estimates = n_raw, kept = nrow(effq), genera_pooled = nrow(p), fdr05 = sum(p$p_adjust < 0.05, na.rm = TRUE),
        up = sum(p$p_adjust < 0.05 & p$estimate > 0, na.rm = TRUE), down = sum(p$p_adjust < 0.05 & p$estimate < 0, na.rm = TRUE),
        pos_controls_fdr05 = sum(pcs$p_adjust < 0.05 & pcs$estimate > 0, na.rm = TRUE),
        median_I2 = round(median(p$I2, na.rm = TRUE), 2), max_weight_over_half = sum(p$max_weight > 0.5, na.rm = TRUE),
        flagged_se = sum(p$n_flagged_se > 0, na.rm = TRUE))
    }
  }
}
S2 <- do.call(rbind, summ); rownames(S2) <- NULL
f <- file.path(odir, "summary.csv"); if (file.exists(f)) { old <- read.csv(f); old <- old[!old$engine %in% S2$engine, ]; S2 <- rbind(old, S2) }
write.csv(S2, f, row.names = FALSE); print(S2)
