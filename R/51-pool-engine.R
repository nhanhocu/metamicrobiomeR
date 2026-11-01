# The lean random-effects engine behind pooling

.rma_fast <- function(y, v, X = NULL, method = c("REML", "DL", "PM"), test = c("z", "knha", "adhoc"), level = 0.95) {
  method <- match.arg(method); test <- match.arg(test)
  k <- length(y)
  scalar <- is.null(X)                                              # intercept only: plain arithmetic, no matrices
  if (scalar) X <- matrix(1, k, 1)
  p <- ncol(X)
  if (k == 1 && p == 1 && is.finite(y) && is.finite(v) && v > 0) {          # one study: its own estimate, no heterogeneity
    se <- sqrt(v); z <- stats::qnorm(1 - (1 - level) / 2)
    return(list(beta = y, se = se, ci_lb = y - z * se, ci_ub = y + z * se, stat = y / se, pval = 2 * stats::pnorm(-abs(y / se)),
                tau2 = 0, I2 = 0, QE = 0, QEp = 1, k = 1L, p = 1L, df = Inf, vcov = matrix(v)))
  }
  if (k <= p || any(!is.finite(y)) || any(!is.finite(v)) || any(v <= 0)) return(NULL)
  if (scalar) return(.rma_scalar(y, v, method, test, level))
  w0 <- 1 / v
  A0inv <- tryCatch(solve(crossprod(X, X * w0)), error = function(e) NULL)
  if (is.null(A0inv)) return(NULL)
  b0 <- A0inv %*% crossprod(X, w0 * y); r0 <- as.numeric(y - X %*% b0)
  QE <- sum(w0 * r0^2)
  trP0 <- sum(w0) - sum(diag(A0inv %*% crossprod(X, X * w0^2)))
  # weighted least squares at a given tau2
  wls <- function(tau2) {
    w <- 1 / (v + tau2); Ainv <- solve(crossprod(X, X * w)); b <- Ainv %*% crossprod(X, w * y)
    list(w = w, Ainv = Ainv, b = b, r = as.numeric(y - X %*% b))
  }
  root <- function(f, f0) .rma_root(f, f0, y, v)
  tau2 <- switch(method,
    DL = max(0, (QE - (k - p)) / trP0),
    PM = { f <- function(t) { s <- wls(t); sum(s$w * s$r^2) - (k - p) }; root(f, QE - (k - p)) },
    REML = {
      f <- function(t) { s <- wls(t); -0.5 * (sum(s$w) - sum(diag(s$Ainv %*% crossprod(X, X * s$w^2)))) + 0.5 * sum(s$w^2 * s$r^2) }
      root(f, f(0))
    })
  s <- wls(tau2)
  q <- sum(s$w * s$r^2) / (k - p)
  scale <- switch(test, z = 1, knha = q, adhoc = max(1, q))
  V <- s$Ainv * scale
  beta <- as.numeric(s$b); se <- sqrt(diag(V)); stat <- beta / se
  df <- if (test == "z") Inf else k - p
  pval <- 2 * stats::pt(-abs(stat), df); crit <- stats::qt(1 - (1 - level) / 2, df)
  vt <- (k - p) / trP0
  out <- list(beta = beta, se = se, ci_lb = beta - crit * se, ci_ub = beta + crit * se, stat = stat, pval = pval,
              tau2 = tau2, I2 = 100 * tau2 / (vt + tau2), QE = QE, QEp = stats::pchisq(QE, k - p, lower.tail = FALSE),
              k = k, p = p, df = df, vcov = V)
  if (p > 1) {
    m <- p - 1; Vm <- V[-1, -1, drop = FALSE]; bm <- beta[-1]
    QM <- as.numeric(crossprod(bm, solve(Vm, bm)))
    out$QM <- if (test == "z") QM else QM / m
    out$QMp <- if (test == "z") stats::pchisq(QM, m, lower.tail = FALSE) else stats::pf(QM / m, m, k - p, lower.tail = FALSE)
  }
  if (p == 1 && k >= 3) {                          # prediction interval of the pooled effect
    cp <- stats::qt(1 - (1 - level) / 2, if (test == "z") Inf else k - p); half <- cp * sqrt(se[1]^2 + tau2)
    out$pi_lb <- beta[1] - half; out$pi_ub <- beta[1] + half
  }
  out
}

# f is decreasing in tau2; its root is the estimate, or 0 if f(0) <= 0
.rma_root <- function(f, f0, y, v) {
  if (f0 <= 0) return(0)
  hi <- max(1, 10 * max(stats::var(y), v)); n <- 0
  while (f(hi) > 0 && n < 12) { hi <- hi * 10; n <- n + 1 }
  if (f(hi) > 0) return(hi)
  stats::uniroot(f, c(0, hi), tol = 1e-10)$root
}

# The intercept-only case with plain arithmetic: same results as the matrix version, several times faster
.rma_scalar <- function(y, v, method, test, level) {
  k <- length(y); w0 <- 1 / v; sw0 <- sum(w0); mu0 <- sum(w0 * y) / sw0
  QE <- sum(w0 * (y - mu0)^2); trP0 <- sw0 - sum(w0^2) / sw0
  state <- function(t) { w <- 1 / (v + t); sw <- sum(w); mu <- sum(w * y) / sw; list(w = w, sw = sw, mu = mu, r = y - mu) }
  tau2 <- switch(method,
    DL = max(0, (QE - (k - 1)) / trP0),
    PM = { f <- function(t) { s <- state(t); sum(s$w * s$r^2) - (k - 1) }; .rma_root(f, QE - (k - 1), y, v) },
    REML = { f <- function(t) { s <- state(t); -0.5 * (s$sw - sum(s$w^2) / s$sw) + 0.5 * sum(s$w^2 * s$r^2) }; .rma_root(f, f(0), y, v) })
  s <- state(tau2)
  q <- sum(s$w * s$r^2) / (k - 1)
  var_b <- switch(test, z = 1, knha = q, adhoc = max(1, q)) / s$sw
  se <- sqrt(var_b); stat <- s$mu / se
  df <- if (test == "z") Inf else k - 1
  crit <- stats::qt(1 - (1 - level) / 2, df); vt <- (k - 1) / trP0
  out <- list(beta = s$mu, se = se, ci_lb = s$mu - crit * se, ci_ub = s$mu + crit * se, stat = stat, pval = 2 * stats::pt(-abs(stat), df),
              tau2 = tau2, I2 = 100 * tau2 / (vt + tau2), QE = QE, QEp = stats::pchisq(QE, k - 1, lower.tail = FALSE), k = k, p = 1L, df = df, vcov = matrix(var_b))
  if (k >= 3) { half <- crit * sqrt(se^2 + tau2); out$pi_lb <- s$mu - half; out$pi_ub <- s$mu + half }
  out
}
