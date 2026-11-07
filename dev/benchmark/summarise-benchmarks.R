o <- commandArgs(TRUE)[1]
rd <- function(p, pat) do.call(rbind, lapply(list.files(p, pattern = pat, full.names = TRUE), read.csv))
x <- rd("bench/long", "^B-.*csv$")
L <- do.call(rbind, lapply(split(x, list(x$scenario, x$engine), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], engine = d$engine[1], fitted = sum(!is.na(d$p)), reject_p05 = round(mean(d$p < 0.05, na.rm = TRUE), 3),
             converged = round(mean(d$converged, na.rm = TRUE), 3), mean_estimate = round(mean(d$estimate, na.rm = TRUE), 3), seconds = round(mean(d$secs), 2))))
write.csv(L[order(L$scenario, L$engine), ], file.path(o, "q2-longitudinal-summary.csv"), row.names = FALSE)
x <- rd("bench/q3", "^Q3-.*csv$")
Q <- do.call(rbind, lapply(split(x, list(x$k, x$tau, x$effect, x$method), drop = TRUE), function(d) { n <- d[!d$diff, ]; t <- d[d$diff, ]
  data.frame(k = d$k[1], tau = d$tau[1], effect = d$effect[1], method = d$method[1], nondiff_reject = round(mean(n$p < 0.05, na.rm = TRUE), 3),
             power = if (nrow(t)) round(mean(t$p < 0.05, na.rm = TRUE), 3) else NA) }))
write.csv(Q[order(Q$effect, Q$tau, Q$k, Q$method), ], file.path(o, "q3-pooling-summary.csv"), row.names = FALSE)
x <- read.csv("bench/q4/Q4-crossplatform.csv")
A <- aggregate(cbind(n_analysed, auc, dir, power, false_pos) ~ engine + design, x, function(v) round(mean(v, na.rm = TRUE), 3), na.action = na.pass)
write.csv(A, file.path(o, "q4-crossplatform-summary.csv"), row.names = FALSE)
file.copy("reana/summary.csv", file.path(o, "reanalysis-summary.csv"), overwrite = TRUE)
