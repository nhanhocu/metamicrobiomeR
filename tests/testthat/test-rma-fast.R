# The lean core must agree with metafor, which stays the reference
# metafor's iterative estimators (REML, PM) stop at a loose tolerance by default, so the reference is run with a tight one:
# that tests the implementation against the exact solution rather than against metafor's rounding.
tight <- list(threshold = 1e-12, maxiter = 5000)
# DL and REML agree with metafor to ~1e-10. metafor solves the Paule-Mandel equation less precisely (its root finder keeps a
# loose tolerance even when asked for a tight one), so PM is compared loosely here and verified against its defining equation below.
tol <- function(m) if (m == "PM") 1e-3 else 1e-6
ref <- function(y, v, method, test, mods = NULL) {
  ctl <- if (method == "DL") list() else tight
  suppressWarnings(if (is.null(mods)) metafor::rma.uni(yi = y, vi = v, method = method, test = test, control = ctl)
                   else metafor::rma.uni(yi = y, vi = v, mods = mods, method = method, test = test, control = ctl))
}
set.seed(123)
cases <- function() {
  out <- list()
  for (k in c(2, 3, 4, 6, 10, 25)) for (het in c(0, 0.05, 0.4, 2)) for (rep in 1:3) {
    v <- runif(k, 0.01, 0.4); y <- rnorm(k, 0.3, sqrt(v + het)); out[[length(out) + 1]] <- list(y = y, v = v, k = k, het = het)
  }
  out
}
cs <- cases()

test_that("the fast core reproduces metafor for the intercept-only model (DL, REML, PM; z, knha, adhoc)", {
  skip_if_not_installed("metafor")
  worst <- 0
  for (cc in cs) for (m in c("DL", "REML", "PM")) for (te in c("z", "knha", "adhoc")) {
    r <- ref(cc$y, cc$v, m, te); f <- metamicrobiomeR:::.rma_fast(cc$y, cc$v, method = m, test = te)
    got <- c(f$beta, f$se, f$ci_lb, f$ci_ub, f$pval, f$tau2, f$QE, f$QEp)
    exp <- c(as.numeric(r$beta), r$se, r$ci.lb, r$ci.ub, r$pval, r$tau2, r$QE, r$QEp)
    d <- max(abs(got - exp) / pmax(1, abs(exp)))
    worst <- max(worst, d)
    if (d > tol(m)) fail(sprintf("k=%d het=%g %s/%s differs by %.2g", cc$k, cc$het, m, te, d))
    # I2 is a percentage (points)
    if (abs(f$I2 - r$I2) > 100 * tol(m)) fail(sprintf("k=%d het=%g %s/%s: I2 %.5f vs %.5f", cc$k, cc$het, m, te, f$I2, r$I2))
  }
  expect_lt(worst, 1e-3)
})

test_that("the prediction interval matches metafor", {
  skip_if_not_installed("metafor")
  n <- 0
  for (cc in cs) if (cc$k >= 3) for (te in c("z", "adhoc")) {
    r <- ref(cc$y, cc$v, "REML", te); pr <- stats::predict(r); f <- metamicrobiomeR:::.rma_fast(cc$y, cc$v, method = "REML", test = te)
    expect_equal(c(f$pi_lb, f$pi_ub), c(pr$pi.lb, pr$pi.ub), tolerance = 1e-6); n <- n + 1
  }
  expect_gt(n, 50)
  expect_null(metamicrobiomeR:::.rma_fast(c(1, 2), c(.1, .1))$pi_lb)                      # not reported with fewer than 3 studies
})

test_that("the fast core reproduces metafor for meta-regression (one and two moderators)", {
  skip_if_not_installed("metafor")
  set.seed(5); worst <- 0
  for (k in c(5, 8, 12)) for (rep in 1:6) for (m in c("DL", "REML", "PM")) for (te in c("z", "knha", "adhoc")) {
    v <- runif(k, 0.02, 0.3); x1 <- rep(c(0, 1), length.out = k); x2 <- rnorm(k)
    y <- rnorm(k, 0.2 + 0.5 * x1, sqrt(v + 0.1))
    for (X in list(cbind(1, x1), cbind(1, x1, x2))) {
      r <- ref(y, v, m, te, mods = X[, -1, drop = FALSE]); f <- metamicrobiomeR:::.rma_fast(y, v, X, method = m, test = te)
      got <- c(f$beta, f$se, f$pval, f$tau2, f$QM, f$QMp); exp <- c(as.numeric(r$beta), r$se, r$pval, r$tau2, r$QM, r$QMp)
      d <- max(abs(got - exp) / pmax(1, abs(exp))); worst <- max(worst, d)
      if (d > tol(m)) fail(sprintf("k=%d p=%d %s/%s differs by %.2g", k, ncol(X), m, te, d))
      if (abs(f$I2 - r$I2) > 100 * tol(m)) fail(sprintf("k=%d p=%d %s/%s: I2 %.5f vs %.5f", k, ncol(X), m, te, f$I2, r$I2))
    }
  }
  expect_lt(worst, 1e-3)
})

test_that("the fast core handles the boundary and refuses what cannot be fitted", {
  r <- metamicrobiomeR:::.rma_fast(c(0.5, 0.5, 0.5), c(.1, .1, .1), method = "REML")
  expect_equal(r$tau2, 0); expect_equal(r$beta, 0.5)
  one <- metamicrobiomeR:::.rma_fast(0.4, 0.04); expect_equal(c(one$beta, one$se, one$tau2, one$I2), c(0.4, 0.2, 0, 0))   # one study: its own estimate
  expect_equal(one$pval, 2 * stats::pnorm(-2))
  expect_null(metamicrobiomeR:::.rma_fast(c(1, 2), c(.1, .1), cbind(1, c(0, 1))))          # not more studies than coefficients
  expect_null(metamicrobiomeR:::.rma_fast(c(1, 2, 3), c(.1, 0, .1)))                        # a zero variance
  expect_null(metamicrobiomeR:::.rma_fast(c(1, NA, 3), c(.1, .1, .1)))
})

test_that("against metafor's default settings the fast core agrees to about four digits", {
  skip_if_not_installed("metafor")
  worst <- 0
  for (cc in cs) for (m in c("DL", "REML", "PM")) {
    r <- suppressWarnings(metafor::rma.uni(yi = cc$y, vi = cc$v, method = m)); f <- metamicrobiomeR:::.rma_fast(cc$y, cc$v, method = m, test = "z")
    worst <- max(worst, abs(f$beta - as.numeric(r$beta)) / max(1, abs(r$beta)), abs(f$se - r$se) / max(1, r$se), abs(f$pval - r$pval))
  }
  expect_lt(worst, 1e-3)
})

test_that("the Paule-Mandel solution satisfies the Paule-Mandel equation exactly", {
  set.seed(7); n <- 0
  for (i in 1:80) {
    k <- sample(3:12, 1); v <- runif(k, .01, .4); y <- rnorm(k, .3, sqrt(v + .2)); f <- metamicrobiomeR:::.rma_fast(y, v, method = "PM")
    if (f$tau2 > 0) { w <- 1 / (v + f$tau2); mu <- sum(w * y) / sum(w); expect_lt(abs(sum(w * (y - mu)^2) - (k - 1)), 1e-7); n <- n + 1 }
  }
  expect_gt(n, 20)
})

test_that("the REML solution satisfies the REML score equation exactly", {
  set.seed(8); n <- 0
  for (i in 1:80) {
    k <- sample(3:12, 1); v <- runif(k, .01, .4); y <- rnorm(k, .3, sqrt(v + .2)); f <- metamicrobiomeR:::.rma_fast(y, v, method = "REML")
    if (f$tau2 > 0) {
      w <- 1 / (v + f$tau2); A <- sum(w); r <- y - sum(w * y) / A
      expect_lt(abs(-0.5 * (sum(w) - sum(w^2) / A) + 0.5 * sum(w^2 * r^2)), 1e-7); n <- n + 1
    }
  }
  expect_gt(n, 20)
})

test_that("the scalar pooling path equals the matrix path", {
  set.seed(77); worst <- 0
  for (k in c(2, 3, 4, 7, 15, 40)) for (rep in 1:8) for (m in c("DL", "REML", "PM")) for (te in c("z", "knha", "adhoc")) {
    v <- runif(k, 0.01, 0.4); y <- rnorm(k, 0.3, sqrt(v + c(0, 0.1, 1)[1 + rep %% 3]))
    a <- metamicrobiomeR:::.rma_fast(y, v, method = m, test = te)                                    # scalar path
    b <- metamicrobiomeR:::.rma_fast(y, v, X = matrix(1, k, 1), method = m, test = te)               # matrix path
    for (nm in c("beta", "se", "ci_lb", "ci_ub", "stat", "pval", "tau2", "I2", "QE", "QEp", "pi_lb", "pi_ub")) {
      x <- a[[nm]]; z <- if (nm %in% c("pi_lb", "pi_ub")) NULL else b[[nm]]
      if (!is.null(z)) worst <- max(worst, abs(x - z) / max(1, abs(z)))
    }
    expect_equal(is.null(a$pi_lb), k < 3)
  }
  expect_lt(worst, 1e-9)
})
