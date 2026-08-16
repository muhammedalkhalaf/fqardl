## ============================================================================
## ADVERSARIAL AUDIT of fqardl 1.0.3
## Every claim published in NEWS.md, the README notice, the GitHub issue and the
## two emails is re-derived here from an independent route. The point is to try
## to BREAK the release, not to confirm it.
## ============================================================================

OUT <- commandArgs(TRUE)[1]
if (is.na(OUT)) OUT <- "adversarial-report.txt"
sink(OUT, split = TRUE)

LIB103 <- .libPaths()[1]
LIB102 <- file.path(tempdir(), "lib102")

PASS <- 0L; FAIL <- 0L; NOTES <- character(0)
ok <- function(lab, cond, detail = "") {
  good <- isTRUE(cond)
  if (good) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  cat(sprintf("  [%s] %s%s\n", if (good) "PASS" else "**FAIL**", lab,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
}
note <- function(x) { NOTES <<- c(NOTES, x); cat("  [note] ", x, "\n", sep = "") }
hr <- function(s) cat("\n", strrep("=", 78), "\n ", s, "\n", strrep("=", 78), "\n", sep = "")

suppressMessages(library(fqardl))
suppressMessages(library(quantreg))
data(macro_data,   package = "fqardl")
data(oil_gdp_data, package = "fqardl")

cat("fqardl version under test:", as.character(packageVersion("fqardl")), "\n")
cat("R:", R.version.string, "\n")

## ===========================================================================
hr("A. PROVENANCE: is the installed package the code that is on GitHub main?")
## ===========================================================================
# The installed source is compared file by file against a fresh clone of main.
CLONE <- file.path(tempdir(), "verify-clone")
unlink(CLONE, recursive = TRUE, force = TRUE)
gitok <- suppressWarnings(system2("git", c("clone", "--depth", "1",
  "https://github.com/muhammedalkhalaf/fqardl.git", shQuote(CLONE)),
  stdout = TRUE, stderr = TRUE))
if (dir.exists(file.path(CLONE, "R"))) {
  local_R <- sort(list.files(file.path(CLONE, "R"), full.names = TRUE, pattern = "[.]R$"))
  inst_R  <- file.path(system.file(package = "fqardl"), "R")
  ok("main contains an R/ directory with 12 sources", length(local_R) == 12,
     paste(length(local_R), "files"))
  ok("main contains NAMESPACE, DESCRIPTION, data, man, tests",
     all(file.exists(file.path(CLONE, c("NAMESPACE", "DESCRIPTION")))) &&
     all(dir.exists(file.path(CLONE, c("data", "man", "tests")))))
  desc <- read.dcf(file.path(CLONE, "DESCRIPTION"))
  ok("DESCRIPTION on main says Version 1.0.3", desc[1, "Version"] == "1.0.3",
     desc[1, "Version"])
  # the decisive line
  src <- readLines(file.path(CLONE, "R", "qardl.R"), warn = FALSE)
  lhs <- grep("y_adj *<-", src, value = TRUE)
  ok("the LHS line on main is the differenced one",
     any(grepl("y\\[max_lag:\\(n - 1\\)\\]", lhs)), trimws(lhs[1]))
  ok("no bare levels LHS remains anywhere on main",
     !any(vapply(list.files(file.path(CLONE, "R"), full.names = TRUE),
                 function(f) any(grepl("^\\s*y_adj <- y\\[\\(max_lag ?\\+ ?1\\):n\\]\\s*$",
                                       readLines(f, warn = FALSE))), logical(1))))
  ok("no mean(t^2) construction remains in code (comments aside)",
     !any(vapply(list.files(file.path(CLONE, "R"), full.names = TRUE),
                 function(f) any(grepl("^\\s*F_stat *<- *mean\\(", readLines(f, warn = FALSE))),
                 logical(1))))
} else {
  note("could not clone the repository; provenance section skipped")
}

## ===========================================================================
hr("B. EXTERNAL ORACLE: rebuild the model by hand and compare")
## ===========================================================================
# fqardl is not allowed to check itself. The same conditional ECM is built from
# scratch with quantreg and the coefficients must agree to machine precision.
y  <- macro_data$gdp
X  <- as.matrix(macro_data[, c("inflation", "interest_rate")])
n  <- length(y); k <- ncol(X)
K  <- 1; p <- 1; q <- 1; m <- max(p, q)
tt <- (m + 1):n
oracle <- data.frame(
  dy      = y[tt] - y[tt - 1],
  y_lag1  = y[tt - 1],
  x1_lag1 = X[tt - 1, 1],
  x2_lag1 = X[tt - 1, 2],
  dx1     = X[tt, 1] - X[tt - 1, 1],
  dx2     = X[tt, 2] - X[tt - 1, 2],
  sin_1   = sin(2 * pi * K * tt / n),
  cos_1   = cos(2 * pi * K * tt / n))
fitO <- rq(dy ~ y_lag1 + x1_lag1 + x2_lag1 + dx1 + dx2 + sin_1 + cos_1,
           data = oracle, tau = 0.5)

set.seed(123)
fitP <- fqardl(gdp ~ inflation + interest_rate, data = macro_data,
               tau = 0.5, max_k = 1, max_p = 1, max_q = 1, verbose = FALSE)
cp <- fitP$qardl_results[[1]]$coefficients
co <- coef(fitO)
cmp <- c("y_lag1", "x1_lag1", "x2_lag1", "dx1", "dx2")
ok("every coefficient matches an independently built ECM",
   isTRUE(all.equal(unname(cp[cmp]), unname(co[cmp]), tolerance = 1e-8)),
   paste(sprintf("%s: pkg %.6f vs oracle %.6f", cmp, cp[cmp], co[cmp]), collapse = "; "))

ok("the long-run multiplier equals -theta/rho computed by hand",
   isTRUE(all.equal(unname(fitP$long_run[1, "inflation"]),
                    unname(-co["x1_lag1"] / co["y_lag1"]), tolerance = 1e-8)))

## ===========================================================================
hr("C. THE WALD STATISTIC: three independent derivations must agree")
## ===========================================================================
r  <- fitP$qardl_results[[1]]
lev <- c("y_lag1", "x1_lag1", "x2_lag1")
pos <- match(lev, names(r$coefficients))
R1  <- matrix(0, length(pos), length(r$coefficients)); R1[cbind(seq_along(pos), pos)] <- 1
Rb  <- R1 %*% r$coefficients
F_manual <- as.numeric(t(Rb) %*% solve(R1 %*% r$vcov %*% t(R1)) %*% Rb) / length(pos)
ok("package F equals the hand-built Wald F",
   isTRUE(all.equal(fitP$bounds_test$F_stat, F_manual, tolerance = 1e-10)),
   sprintf("pkg %.6f vs manual %.6f", fitP$bounds_test$F_stat, F_manual))
# third route: quadratic form written out elementwise
V  <- r$vcov[pos, pos]
b  <- r$coefficients[pos]
F_qf <- as.numeric(crossprod(b, solve(V, b))) / length(pos)
ok("and equals the elementwise quadratic form b'V^{-1}b/m",
   isTRUE(all.equal(fitP$bounds_test$F_stat, F_qf, tolerance = 1e-10)))
# it must NOT equal the discarded construction
ok("and is NOT mean(t^2)",
   !isTRUE(all.equal(fitP$bounds_test$F_stat, mean(r$t_statistics[lev]^2), tolerance = 1e-6)),
   sprintf("Wald %.3f vs mean(t^2) %.3f", fitP$bounds_test$F_stat,
           mean(r$t_statistics[lev]^2)))
# the divisor must be the number of restrictions, k+1
ok("the divisor is k+1", fitP$bounds_test$n_restrictions == k + 1)

## ===========================================================================
hr("D. SIZE UNDER THE NULL: does the bounds test still over-reject?")
## ===========================================================================
# Two INDEPENDENT random walks. There is no cointegration, so a correctly sized
# 5 percent test should exceed the I(1) upper bound about 5 percent of the time.
# Version 1.0.2 exceeded it essentially always. This is the sharpest test of the
# whole repair.
Rn <- 150
rej <- rej_t <- rep(NA, Rn)
for (i in seq_len(Rn)) {
  set.seed(50000 + i)
  nn <- 150
  d <- data.frame(y = cumsum(rnorm(nn)), x1 = cumsum(rnorm(nn)))
  set.seed(3)
  f <- try(fqardl(y ~ x1, data = d, tau = 0.5, max_k = 1, max_p = 2, max_q = 2,
                  verbose = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) next
  cv <- fqardl:::get_pss_critical_values(1, case = 3)
  rej[i]   <- f$bounds_test$F_stat > cv$F_upper[2]     # 5 percent I(1) bound
  rej_t[i] <- f$bounds_test$t_stat < cv$t_upper[2]
}
rr  <- mean(rej,   na.rm = TRUE)
rrt <- mean(rej_t, na.rm = TRUE)
cat(sprintf("  usable replications: %d\n", sum(!is.na(rej))))
cat(sprintf("  F exceeds the 5%% I(1) bound in %.1f%% of samples\n", 100 * rr))
cat(sprintf("  t falls below the 5%% I(1) bound in %.1f%% of samples\n", 100 * rrt))
ok("empirical size of the F route is not grossly inflated (< 20%)", rr < 0.20,
   sprintf("%.1f%%", 100 * rr))
ok("empirical size of the t route is not grossly inflated (< 20%)", rrt < 0.20,
   sprintf("%.1f%%", 100 * rrt))
if (rr > 0.10) note(sprintf(
  "size is %.1f%% against a nominal 5%%: the PSS bounds are derived for the conditional mean and for OLS, not for a bootstrapped quantile covariance, so some over-rejection is expected and is already documented as a limitation.", 100 * rr))

## ===========================================================================
hr("E. RECOVERY ACROSS THE PARAMETER SPACE, not just one design")
## ===========================================================================
grid <- expand.grid(rho = c(-0.10, -0.35, -0.70), beta = c(-1.5, 0.8, 2.0))
res <- NULL
for (g in seq_len(nrow(grid))) {
  rho <- grid$rho[g]; bet <- grid$beta[g]
  est <- rep(NA_real_, 12); ect <- rep(NA_real_, 12)
  for (i in 1:12) {
    set.seed(9000 + 100 * g + i); nn <- 300
    x <- cumsum(rnorm(nn)); yy <- numeric(nn); yy[1] <- bet * x[1]
    for (t in 2:nn)
      yy[t] <- yy[t-1] + rho * (yy[t-1] - 2 - bet * x[t-1]) +
               0.8 * (x[t] - x[t-1]) + rnorm(1, 0, .5)
    set.seed(3)
    f <- try(fqardl(y ~ x1, data = data.frame(y = yy, x1 = x), tau = 0.5,
                    max_k = 1, max_p = 2, max_q = 2, verbose = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) next
    est[i] <- f$long_run[1, "x1"]
    ect[i] <- f$qardl_results[[1]]$coefficients["y_lag1"]
  }
  res <- rbind(res, data.frame(rho_true = rho, beta_true = bet,
                               ect_med = median(ect, na.rm = TRUE),
                               beta_med = median(est, na.rm = TRUE),
                               sign_ok = mean(sign(est) == sign(bet), na.rm = TRUE)))
}
print(round(res, 3), row.names = FALSE)
ok("the long-run sign is right in every cell of the grid", all(res$sign_ok == 1))
ok("the long-run estimate is within 0.15 of truth in every cell",
   all(abs(res$beta_med - res$beta_true) < 0.15),
   sprintf("max abs error %.3f", max(abs(res$beta_med - res$beta_true))))
ok("the ECT is negative in every cell", all(res$ect_med < 0))
ok("the ECT tracks the true adjustment speed (correlation > 0.9)",
   cor(res$ect_med, res$rho_true) > 0.9,
   sprintf("cor = %.3f", cor(res$ect_med, res$rho_true)))

## ===========================================================================
hr("F. DIMENSION STRESS: k != q, several p, and the level-term selection")
## ===========================================================================
# The qardlr sibling had a lag-major / variable-major mismatch that only bites
# when k differs from q. Check that fqardl selects exactly y_lag1 plus every
# x*_lag1, whatever k, p and q are.
bad <- character(0)
for (kk in 1:4) for (pp in 1:3) for (qq in 1:3) {
  set.seed(400 + 10 * kk + pp + qq); nn <- 200
  Xg <- matrix(cumsum(rnorm(nn * kk)), nn, kk)
  colnames(Xg) <- paste0("v", seq_len(kk))
  yg <- numeric(nn)
  for (t in 2:nn) yg[t] <- yg[t-1] - 0.4 * (yg[t-1] - 0.7 * Xg[t-1, 1]) + rnorm(1)
  dd <- data.frame(y = yg, Xg)
  fo <- as.formula(paste("y ~", paste(colnames(Xg), collapse = "+")))
  set.seed(3)
  f <- try(fqardl(fo, data = dd, tau = 0.5, max_k = 1,
                  max_p = pp, max_q = qq, verbose = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) { bad <- c(bad, sprintf("k=%d p=%d q=%d ERROR", kk, pp, qq)); next }
  nm  <- names(f$qardl_results[[1]]$coefficients)
  want <- c("y_lag1", paste0("x", seq_len(kk), "_lag1"))
  if (!all(want %in% nm)) bad <- c(bad, sprintf("k=%d p=%d q=%d missing level terms", kk, pp, qq))
  if (f$bounds_test$n_restrictions != kk + 1)
    bad <- c(bad, sprintf("k=%d p=%d q=%d restrictions=%d expected %d", kk, pp, qq,
                          f$bounds_test$n_restrictions, kk + 1))
  if (!is.finite(f$bounds_test$F_stat))
    bad <- c(bad, sprintf("k=%d p=%d q=%d F not finite", kk, pp, qq))
  if (!is.finite(f$bounds_test$t_stat))
    bad <- c(bad, sprintf("k=%d p=%d q=%d t not finite", kk, pp, qq))
  # long-run must equal -theta/rho for EVERY covariate, not just the first
  cf <- f$qardl_results[[1]]$coefficients
  lr_hand <- -cf[paste0("x", seq_len(kk), "_lag1")] / cf["y_lag1"]
  if (!isTRUE(all.equal(unname(f$long_run[1, ]), unname(lr_hand), tolerance = 1e-8)))
    bad <- c(bad, sprintf("k=%d p=%d q=%d long-run mismatch", kk, pp, qq))
}
ok("36 combinations of k, p and q all behave", length(bad) == 0,
   if (length(bad)) paste(head(bad, 6), collapse = " | ") else "k=1..4, p=1..3, q=1..3")

## ===========================================================================
hr("G. FALSIFICATION: undo the fix and the old symptom must return")
## ===========================================================================
# A repair you cannot un-do is a repair you have not localised.
ORIG <- fqardl:::build_ardl_design
BROKEN <- function(y, X, fourier, p, q) {
  d <- ORIG(y, X, fourier, p, q); if (is.null(d)) return(NULL)
  mm <- d$max_lag; N <- length(y); d$y <- y[(mm + 1):N]; d      # the 1.0.2 line
}
assignInNamespace("build_ardl_design", BROKEN, ns = "fqardl")
set.seed(123)
fB <- fqardl(gdp ~ inflation + interest_rate, data = macro_data,
             tau = c(.10,.25,.50,.75,.90), max_k = 2, verbose = FALSE)
assignInNamespace("build_ardl_design", ORIG, ns = "fqardl")
set.seed(123)
fG <- fqardl(gdp ~ inflation + interest_rate, data = macro_data,
             tau = c(.10,.25,.50,.75,.90), max_k = 2, verbose = FALSE)
cB <- vapply(fB$qardl_results, function(z) unname(z$coefficients["y_lag1"]), numeric(1))
cG <- vapply(fG$qardl_results, function(z) unname(z$coefficients["y_lag1"]), numeric(1))
print(round(data.frame(tau = fG$tau, broken = cB, fixed = cG, diff = cB - cG), 4),
      row.names = FALSE)
ok("re-introducing the old line makes every coefficient positive again", all(cB > 0))
ok("the difference is exactly 1 at every quantile (phi = 1 + rho)",
   isTRUE(all.equal(unname(cB - cG), rep(1, length(cB)), tolerance = 1e-9)))
# NOTE: macro_data is not cointegrated (the corrected F is 2.34, below the I(0)
# bound), so a slightly positive ECT at an extreme quantile is a legitimate
# outcome, not a defect. Asserting all(cG < 0) here was an error in this audit
# script, corrected below: what must hold is that the fix removes the systematic
# near-unity positive coefficient, not that every quantile turns negative.
ok("the fix removes the near-unity coefficients (all now far below 0.5)", all(cG < 0.5))
ok("on a genuinely cointegrated design every quantile is negative",
   { set.seed(77); nn <- 300
     xx <- cumsum(rnorm(nn)); yy <- numeric(nn); yy[1] <- 0.8*xx[1]
     for (t in 2:nn) yy[t] <- yy[t-1] - 0.45*(yy[t-1] - 0.8*xx[t-1]) +
                              0.8*(xx[t]-xx[t-1]) + rnorm(1,0,.5)
     set.seed(3)
     fc <- fqardl(y ~ x1, data = data.frame(y=yy, x1=xx), tau = c(.1,.25,.5,.75,.9),
                  max_k = 1, max_p = 2, max_q = 2, verbose = FALSE)
     all(vapply(fc$qardl_results, function(z) z$coefficients["y_lag1"], numeric(1)) < 0) })

## ===========================================================================
hr("H. THE PUBLISHED NUMBERS: each one re-derived")
## ===========================================================================
cat("  claim: macro_data bounds F falls from 65.57 to 2.34\n")
# 1.0.2's F on this design = mean(t^2) over the level terms of the LEVELS fit
rB  <- fB$qardl_results[["tau_0.5"]]
levB <- c("y_lag1", "x1_lag1", "x2_lag1")
F102 <- mean(rB$t_statistics[levB]^2)
F103 <- fG$bounds_test$F_stat
cat(sprintf("    reconstructed 1.0.2 F = %.2f ; 1.0.3 F = %.2f\n", F102, F103))
ok("the 1.0.2 figure reproduces as 65.57 to within rounding", abs(F102 - 65.57) < 0.6,
   sprintf("%.2f", F102))
ok("the 1.0.3 figure reproduces as 2.34 to within rounding", abs(F103 - 2.34) < 0.05,
   sprintf("%.2f", F103))
cvk2 <- fqardl:::get_pss_critical_values(2, case = 3)
ok("2.34 is indeed below the 5 percent I(0) bound of 3.79",
   F103 < cvk2$F_lower[2], sprintf("bound %.2f", cvk2$F_lower[2]))
ok("65.57 is indeed above the 5 percent I(1) bound of 4.85",
   F102 > cvk2$F_upper[2], sprintf("bound %.2f", cvk2$F_upper[2]))
cat("  1.0.2 additionally used a WRONG bound of 4.38 at k = 2; the paper gives 4.85.\n")

cat("\n  claim: fnardl on oil_gdp_data moves from +0.894 / t=+27.08 to -0.106 / t=-3.21, F=4.47\n")
set.seed(11)
fn <- fnardl(gdp ~ oil_price, data = oil_gdp_data, max_k = 1, verbose = FALSE)
cat(sprintf("    fnardl now: ECT %+.4f  t %+.3f  F %.3f\n",
            fn$nardl_result$coefficients["y_lag1"],
            fn$nardl_result$t_statistics["y_lag1"], fn$bounds_test$F_stat))
ok("fnardl ECT reproduces as -0.106", abs(fn$nardl_result$coefficients["y_lag1"] + 0.106) < 0.01)
ok("fnardl t reproduces as -3.21", abs(fn$nardl_result$t_statistics["y_lag1"] + 3.21) < 0.05)
ok("fnardl F reproduces as 4.47", abs(fn$bounds_test$F_stat - 4.47) < 0.05)

cat("\n  claim: mtnardl returns a finite t and a negative ECT\n")
set.seed(11)
mt <- mtnardl(gdp ~ oil_price, data = oil_gdp_data, verbose = FALSE)
cat(sprintf("    mtnardl now: ECT %+.4f  t %+.3f  F %.3f\n",
            mt$model$coefficients["y_lag1"], mt$model$t_statistics["y_lag1"],
            mt$bounds_test$F_stat))
ok("mtnardl ECT negative", unname(mt$model$coefficients["y_lag1"]) < 0)
ok("mtnardl t finite and negative", is.finite(mt$bounds_test$t_stat) && mt$bounds_test$t_stat < 0)

## ===========================================================================
hr("I. THE CRITICAL VALUE TABLE, checked against outside anchors")
## ===========================================================================
# Anchor 1: PSS say the I(0) column of Table CII is constant in k.
tl <- vapply(1:10, function(kk) fqardl:::get_pss_critical_values(kk, 3)$t_lower[2], numeric(1))
ok("the t table's I(0) bound does not vary with k", length(unique(round(tl, 6))) == 1,
   sprintf("%.2f", tl[1]))
# Anchor 2: at k = 0 the Case III t bound equals the Dickey-Fuller value -2.86.
ok("that constant equals the Dickey-Fuller value -2.86", abs(tl[1] + 2.86) < 1e-9)
# Anchor 3: the I(1) bound must fall monotonically in k, and stay below I(0).
tu <- vapply(1:10, function(kk) fqardl:::get_pss_critical_values(kk, 3)$t_upper[2], numeric(1))
ok("the t I(1) bound is monotone decreasing in k", all(diff(tu) < 0))
ok("the t I(1) bound is always below the I(0) bound", all(tu < tl))
# Anchor 4: F bounds must be positive, ordered, and decreasing in k.
Fl <- vapply(1:10, function(kk) fqardl:::get_pss_critical_values(kk, 3)$F_lower[2], numeric(1))
Fu <- vapply(1:10, function(kk) fqardl:::get_pss_critical_values(kk, 3)$F_upper[2], numeric(1))
ok("F I(1) exceeds F I(0) at every k", all(Fu > Fl))
ok("both F bounds decrease in k", all(diff(Fl) < 0) && all(diff(Fu) < 0))
ok("the unverified t entries are NA, not invented",
   all(is.na(vapply(1:10, function(kk)
     fqardl:::get_pss_critical_values(kk, 3)$t_upper[1], numeric(1)))))

## the 1 and 10 percent points must now be absent, not invented
ok("the 1 and 10 percent F bounds are NA, not numbers of unknown provenance",
   all(vapply(1:10, function(kk) {
     cv <- fqardl:::get_pss_critical_values(kk, 3)
     all(is.na(c(cv$F_lower[c(1,3)], cv$F_upper[c(1,3)])))
   }, logical(1))))
ok("the table still varies with k beyond k = 6 (1.0.2 collapsed to one row)",
   fqardl:::get_pss_critical_values(7, 3)$F_upper[2] >
   fqardl:::get_pss_critical_values(10, 3)$F_upper[2])
ok("the shipped 5 percent values differ from the wrong 1.0.2 table at k = 2",
   !isTRUE(all.equal(c(fqardl:::get_pss_critical_values(2,3)$F_lower[2],
                       fqardl:::get_pss_critical_values(2,3)$F_upper[2]),
                     c(3.55, 4.38), tolerance = 1e-8)))

## ===========================================================================
hr("J. THINGS THE REPAIR ITSELF MIGHT HAVE BROKEN")
## ===========================================================================
ok("summary() prints without error",
   !inherits(try(capture.output(summary(fG)), silent = TRUE), "try-error"))
ok("print() prints without error",
   !inherits(try(capture.output(print(fG)), silent = TRUE), "try-error"))
ok("the 1.0.2 return shape survives: $bounds_test$F_stat is still a scalar",
   is.numeric(fG$bounds_test$F_stat) && length(fG$bounds_test$F_stat) == 1)
ok("$bounds_test$decision is still a single string",
   is.character(fG$bounds_test$decision) && length(fG$bounds_test$decision) == 1)
ok("the new per-quantile table is present and has one row per tau",
   nrow(fG$bounds_test$summary) == length(fG$tau))
ok("vcov is returned and is square and symmetric",
   { V <- fG$qardl_results[[1]]$vcov
     !is.null(V) && nrow(V) == ncol(V) && isTRUE(all.equal(V, t(V), tolerance = 1e-8)) })
ok("vcov is positive semi-definite",
   min(eigen(fG$qardl_results[[1]]$vcov, only.values = TRUE)$values) > -1e-8)
w <- tryCatch({ quantile_wald_test(fG$qardl_results, "y_lag1"); "no warning" },
              warning = function(w) conditionMessage(w))
ok("quantile_wald_test still runs and warns about its limitation",
   grepl("independence across quantiles", w))
bb <- tryCatch({ suppressWarnings(bootstrap_bounds_test(
        macro_data$gdp, as.matrix(macro_data[, "inflation", drop = FALSE]),
        generate_fourier_terms(nrow(macro_data), 1), 1, 1, 0.5, 3, n_boot = 20)) },
        error = function(e) e)
ok("bootstrap_bounds_test still runs end to end", !inherits(bb, "error"),
   if (inherits(bb, "error")) conditionMessage(bb) else
     sprintf("p_F = %.3f", bb$p_value_F))
ok("unit root routines are untouched and still run",
   !inherits(try(fourier_adf_test(macro_data$gdp, max_freq = 2, verbose = FALSE), silent = TRUE), "try-error"))
ok("case = 3 still works while other cases error",
   !inherits(try(fqardl(gdp ~ inflation, data = macro_data, tau = .5, max_k = 1,
                        case = 3, verbose = FALSE), silent = TRUE), "try-error") &&
   inherits(try(fqardl(gdp ~ inflation, data = macro_data, tau = .5, max_k = 1,
                        case = 4, verbose = FALSE), silent = TRUE), "try-error"))

## ===========================================================================
hr("K. BOUNDARY AND DEGENERATE INPUTS")
## ===========================================================================
edge <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  cat(sprintf("  %-46s -> %s\n", lab,
      if (inherits(r, "error")) paste("error:", substr(conditionMessage(r), 1, 60))
      else if (inherits(r, "warning")) paste("warning:", substr(conditionMessage(r), 1, 50))
      else "ran"))
  invisible(r)
}
set.seed(5); nn <- 60
dsmall <- data.frame(y = cumsum(rnorm(nn)), x1 = cumsum(rnorm(nn)))
edge("p = 1, q = 1 (the Stata original crashes here)",
     fqardl(y ~ x1, data = dsmall, tau = .5, max_k = 1, max_p = 1, max_q = 1, verbose = FALSE))
edge("tau very low (0.05)",
     fqardl(y ~ x1, data = dsmall, tau = .05, max_k = 1, max_p = 2, max_q = 2, verbose = FALSE))
edge("tau very high (0.95)",
     fqardl(y ~ x1, data = dsmall, tau = .95, max_k = 1, max_p = 2, max_q = 2, verbose = FALSE))
edge("a constant regressor",
     fqardl(y ~ x1, data = data.frame(y = dsmall$y, x1 = rep(1, nn)),
            tau = .5, max_k = 1, verbose = FALSE))
edge("perfectly collinear regressors",
     fqardl(y ~ x1 + x2, data = data.frame(y = dsmall$y, x1 = dsmall$x1,
            x2 = 2 * dsmall$x1), tau = .5, max_k = 1, verbose = FALSE))
edge("missing values present",
     { d2 <- dsmall; d2$x1[c(5, 20)] <- NA
       fqardl(y ~ x1, data = d2, tau = .5, max_k = 1, verbose = FALSE) })
note("Errors above are acceptable outcomes; silent wrong numbers are not. Read each line.")

## ===========================================================================
hr("L. 1.0.2 SIDE BY SIDE, if it can be installed")
## ===========================================================================
if (!dir.exists(LIB102)) dir.create(LIB102, recursive = TRUE)
got <- try(suppressWarnings(install.packages(
  "https://cran.r-project.org/src/contrib/Archive/fqardl/fqardl_1.0.2.tar.gz",
  lib = LIB102, repos = NULL, type = "source", quiet = TRUE)), silent = TRUE)
if (!requireNamespace("fqardl", lib.loc = LIB102, quietly = TRUE)) {
  got <- try(suppressWarnings(install.packages(
    "https://cran.r-project.org/src/contrib/fqardl_1.0.2.tar.gz",
    lib = LIB102, repos = NULL, type = "source", quiet = TRUE)), silent = TRUE)
}
if (dir.exists(file.path(LIB102, "fqardl"))) {
  cat("  fqardl 1.0.2 installed side by side; running it in a separate process\n")
  scr <- file.path(tempdir(), "run102.R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', gsub("\\\\", "/", LIB102)),
    'suppressMessages(library(fqardl))',
    'data(macro_data, package="fqardl")',
    'set.seed(123)',
    'f <- fqardl(gdp ~ inflation + interest_rate, data=macro_data,',
    '            tau=c(.10,.25,.50,.75,.90), max_k=2, verbose=FALSE)',
    'cat("VER", as.character(packageVersion("fqardl")), "\\n")',
    'cat("F", f$bounds_test$F_stat, "\\n")',
    'cat("ECT", sapply(f$qardl_results, function(z) z$coefficients["y_lag1"]), "\\n")',
    'cat("LR", as.vector(f$long_run), "\\n")'), scr)
  o <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(scr),
                                stdout = TRUE, stderr = TRUE))
  cat(paste("   ", o, collapse = "\n"), "\n")
  fline <- grep("^F ", o, value = TRUE)
  if (length(fline)) {
    F102real <- as.numeric(strsplit(fline, " ")[[1]][2])
    ok("the real 1.0.2 F on macro_data matches the 65.57 we published",
       abs(F102real - 65.57) < 0.6, sprintf("%.2f", F102real))
  }
} else {
  note("fqardl 1.0.2 could not be installed here (CRAN archive unreachable); section L relies on the reconstruction in section G instead, which reproduces the same design.")
}

## ===========================================================================
hr("SUMMARY")
## ===========================================================================
cat(sprintf("\n  checks passed : %d\n  checks FAILED : %d\n", PASS, FAIL))
if (length(NOTES)) {
  cat("\n  notes requiring a human decision:\n")
  for (nn2 in NOTES) cat("   - ", nn2, "\n", sep = "")
}
cat(sprintf("\n  VERDICT: %s\n", if (FAIL == 0)
  "no counter-evidence found against the 1.0.3 claims" else
  "AT LEAST ONE CLAIM DID NOT SURVIVE. Do not submit to CRAN until resolved."))
sink()
