# =====================================================================
# Probe: can an mcp (JAGS) fit be driven by rb_loo's engines?
#
# Decides ONE question before any extractor is written: do mcp's
# accessors line up, exactly, with what .rb_engine needs?
#
#   Lf    S x N conditional log-likelihood      <- fit$mcmc_loglik
#   etaF  S x N mean with the RE zeroed         <- fitted(varying = FALSE)
#   U     S x G random-effect draws             <- as_draws (naming discovered)
#   sigma S residual SD draws                   <- as_draws
#   gidx  N group index                         <- data
#
# The checks are alignment checks. Each one, if it passes, retires a
# whole class of silent-wrong-answer bugs (draw order, observation
# order, group mapping, what `varying = FALSE` actually means).
#
# Prereq: JAGS (external) + rjags + mcp(dev).
#   install.packages("rjags")
#   remotes::install_github("lindeloev/mcp@dev")
#
# Prints verdicts, dimensions and max abs differences only. No rows.
# =====================================================================

TOL <- 1e-8
.results <- new.env(parent = emptyenv()); .results$tab <- character(0)
ok  <- function(label, pass, detail = "") {
  .results$tab[label] <- if (isTRUE(pass)) "PASS" else "FAIL"
  cat(sprintf("  [%s] %-52s %s\n", if (isTRUE(pass)) "PASS" else "FAIL",
              label, detail))
}
hdr <- function(x) cat("\n== ", x, " ", strrep("=", max(0, 58 - nchar(x))), "\n", sep = "")

hdr("A. environment")
for (p in c("mcp", "rjags", "coda", "posterior")) {
  have <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-10s %-6s %s\n", p, have,
              if (have) as.character(packageVersion(p)) else ""))
}
if (!requireNamespace("mcp", quietly = TRUE) ||
    !requireNamespace("rjags", quietly = TRUE))
  stop("mcp/rjags not available - install JAGS first, then the R packages.")
suppressMessages({ library(mcp); library(posterior) })
options(mc.cores = 1)

# ---------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------
# mcp's draws: prefer the dev as_draws_*; fall back to the mcmc.list.
get_draws <- function(fit) {
  d <- tryCatch(posterior::as_draws_matrix(fit),
                error = function(e)
                  posterior::as_draws_matrix(fit$mcmc_post))
  as.matrix(d)
}
# The pointwise log-lik lives at fit$loglik and is populated lazily by
# add_loglik(). varying = TRUE gives the CONDITIONAL log-likelihood (random
# effects included), which is the Lf .rb_engine needs. Note that
# loo.mcpfit() also defaults to varying = TRUE, so mcp's shipped loo() is
# the conditional LOO -- the estimand rb_loo replaces.
loglik_matrix <- function(fit) {
  add_ll <- if ("add_loglik" %in% getNamespaceExports("mcp"))
    mcp::add_loglik else mcp:::add_loglik
  if (is.null(fit$loglik)) fit <- add_ll(fit, varying = TRUE)
  ll <- fit$loglik
  if (is.null(ll)) stop("could not populate fit$loglik via add_loglik()")
  ll <- as.matrix(ll)
  # rownames carry the chain id; keep them for the ordering check
  attr(ll, "chain_id") <- suppressWarnings(as.integer(rownames(ll)))
  ll
}
fitted_mat <- function(fit, varying) {
  m <- fitted(fit, summary = FALSE, varying = varying,
              samples_format = "matrix")
  matrix(as.numeric(as.matrix(m)), nrow = nrow(as.matrix(m)))
}

# =====================================================================
# B. IN-SCOPE CASE: varying intercept on the FIRST segment.
#    Z is a column of ones and enters additively, so this is the exact
#    Gaussian random-intercept downdate rb_loo already implements.
# =====================================================================
hdr("B. varying intercept, segment 1 (expected IN scope)")

set.seed(42)
J <- 12; nj <- 8
id  <- rep(seq_len(J), each = nj)
x   <- rep(seq(0, 10, length.out = nj), J)
u   <- rnorm(J, 0, 2)                      # subject intercept deviations
cp  <- 5
mu  <- 10 + u[id] + ifelse(x > cp, -1.2 * (x - cp), 0)
d   <- data.frame(y = mu + rnorm(length(x), 0, 0.7), x = x,
                  id = factor(id))

fit_b <- mcp(list(y ~ 1 + (1 | id), 1 ~ 0 + x), data = d,
             chains = 2, iter = 6000, adapt = 2000)

dm   <- get_draws(fit_b)
vars <- colnames(dm)
cat("  variables:", paste(utils::head(vars, 40), collapse = ", "), "\n")

S <- nrow(dm)
N <- nrow(d)
f_on  <- fitted_mat(fit_b, varying = TRUE)
f_off <- fitted_mat(fit_b, varying = FALSE)
ok("fitted(varying=TRUE) is S x N", all(dim(f_on) == c(S, N)),
   sprintf("got %d x %d, expected %d x %d", nrow(f_on), ncol(f_on), S, N))
ok("fitted(varying=FALSE) matches that shape",
   all(dim(f_off) == dim(f_on)))

# --- discover the RE and sigma parameter names (do not hard-code) ----
re_cols <- grep("^Intercept_1_id\\[", vars, value = TRUE)
if (!length(re_cols)) re_cols <- grep("_id\\[", vars, value = TRUE)
sd_col  <- grep("_id_sd$|^sigma_u$", vars, value = TRUE)
sig_col <- grep("^sigma_1$|^sigma$", vars, value = TRUE)
cat("  RE draws :", length(re_cols), "cols e.g.",
    paste(utils::head(re_cols, 3), collapse = ", "), "\n")
cat("  RE sd    :", paste(sd_col, collapse = ", "),
    "| residual sd:", paste(sig_col, collapse = ", "), "\n")
ok("one RE column per group", length(re_cols) == J,
   sprintf("%d cols vs G = %d", length(re_cols), J))
ok("exactly one RE-SD and one residual-SD parameter",
   length(sd_col) == 1L && length(sig_col) == 1L)

# --- CHECK 1: does varying=FALSE zero the deviations, and do draw ----
# order and the group->observation mapping agree?
# fitted(TRUE) - fitted(FALSE) must equal u_draws[s, group(i)] exactly.
U    <- dm[, re_cols, drop = FALSE]                  # S x G
gidx <- as.integer(d$id)
delta <- f_on - f_off
gap   <- max(abs(delta - U[, gidx, drop = FALSE]))
ok("fitted(TRUE)-fitted(FALSE) == RE draw for that obs", gap < TOL,
   sprintf("max |diff| = %.3e", gap))
if (gap >= TOL) {
  # diagnose: is it a column-order (group label) problem only?
  best <- min(sapply(seq_len(S), function(s)
    max(abs(delta[s, ] - U[s, gidx]))))
  cat("     note: best per-draw gap", sprintf("%.3e", best),
      "- if small, suspect group-label order, not draw order\n")
}

# --- CHECK 2: reconstruct the conditional log-lik from the pieces ----
# This is the decisive one. If dnorm(y, fitted_with_RE, sigma) reproduces
# mcmc_loglik, then observation order, draw order and the sigma name are
# ALL confirmed simultaneously.
Lf <- loglik_matrix(fit_b)
ok("fit$loglik is S x N", all(dim(Lf) == c(S, N)),
   sprintf("got %d x %d", nrow(Lf), ncol(Lf)))
cid <- attr(Lf, "chain_id")
if (!is.null(cid) && !anyNA(cid))
  ok("loglik chains stacked in draw order (not interleaved)",
     !is.unsorted(cid), sprintf("chains: %s", paste(unique(cid), collapse = ",")))
if (all(dim(Lf) == c(S, N)) && length(sig_col) == 1L) {
  sig  <- as.numeric(dm[, sig_col])
  Lhat <- matrix(0, S, N)
  for (i in seq_len(N))
    Lhat[, i] <- dnorm(d$y[i], f_on[, i], sig, log = TRUE)
  dev <- max(abs(Lhat - Lf))
  ok("dnorm(y, fitted, sigma) reproduces mcmc_loglik", dev < 1e-6,
     sprintf("max |diff| = %.3e", dev))
}

# --- CHECK 3: is the RE prior an INDEPENDENT N(0, sd^2)? -------------
# rb_loo's downdate assumes u_j ~ N(0, sigma_u^2) independently. A
# sum-to-zero (or otherwise centred) constraint would invalidate it.
rs   <- rowSums(U)
free <- sd(rs) / (mean(as.numeric(dm[, sd_col])) * sqrt(J))
ok("RE draws are NOT sum-to-zero constrained", free > 0.5,
   sprintf("sd(rowSums)/expected = %.3f (near 0 => constrained)", free))
cat(sprintf("     mean(RE sd draw) = %.3f (simulated 2.0)\n",
            mean(as.numeric(dm[, sd_col]))))

# =====================================================================
# C. OUT-OF-SCOPE CASE: varying CHANGE POINT.
#    Recorded for the bigger decision: naming, and whether the prior is
#    truncated (which a callable-mean quadrature would have to respect).
# =====================================================================
hdr("C. varying change point (expected OUT of scope)")

cp_j <- 5 + rnorm(J, 0, 1.0)
mu2  <- 10 + ifelse(x > cp_j[id], -1.2 * (x - cp_j[id]), 0)
d2   <- data.frame(y = mu2 + rnorm(length(x), 0, 0.7), x = x,
                   id = factor(id))
fit_c <- mcp(list(y ~ 1, 1 + (1 | id) ~ 0 + x), data = d2,
             chains = 2, iter = 6000, adapt = 2000)

dm2   <- get_draws(fit_c)
cp_re <- grep("^cp_1_id\\[", colnames(dm2), value = TRUE)
if (!length(cp_re)) cp_re <- grep("cp.*_id\\[", colnames(dm2), value = TRUE)
cat("  cp RE cols:", length(cp_re), "e.g.",
    paste(utils::head(cp_re, 3), collapse = ", "), "\n")

if (length(cp_re)) {
  CP <- dm2[, cp_re, drop = FALSE]
  rs2  <- rowSums(CP)
  cat(sprintf("     sd(rowSums of cp deviations) = %.4f  (0 => zero-centred\n",
              sd(rs2)))
  cat("     within draw, which breaks the independent-N(0,s^2) prior)\n")
  cat(sprintf("     cp deviation range: [%.3f, %.3f]; x range [%.1f, %.1f]\n",
              min(CP), max(CP), min(d2$x), max(d2$x)))
  # truncation leaves a hard edge: deviations cannot cross the segment bounds
  ok("cp deviations stay inside the x range (truncation active)",
     all(CP + mean(dm2[, "cp_1"]) >= min(d2$x)) &&
     all(CP + mean(dm2[, "cp_1"]) <= max(d2$x)))
}
# does the varying=FALSE accessor still behave for a varying cp?
f2_on  <- tryCatch(fitted_mat(fit_c, varying = TRUE),  error = function(e) NULL)
f2_off <- tryCatch(fitted_mat(fit_c, varying = FALSE), error = function(e) NULL)
ok("fitted() works for a varying-cp model",
   !is.null(f2_on) && !is.null(f2_off))
if (!is.null(f2_on) && !is.null(f2_off)) {
  # Additivity test. If the RE entered as a plain shift, then within one
  # draw and one group every observation would move by the SAME amount,
  # so the within-(draw x group) spread of (on - off) would be 0. A
  # varying change point moves only the observations near/after the
  # break, so the spread is large. (Part B's CHECK 1 is the positive
  # control for this: there the spread is 0 by construction.)
  dlt2  <- f2_on - f2_off
  gidx2 <- as.integer(d2$id)
  spread <- max(vapply(split(seq_len(ncol(dlt2)), gidx2), function(cols)
    max(apply(dlt2[, cols, drop = FALSE], 1, function(v) diff(range(v)))),
    numeric(1)))
  cat(sprintf("     within-(draw x group) spread of the RE shift: %.3e\n",
              spread))
  cat("     (~0 => additive, closed-form downdate applies;\n")
  cat("      large => non-additive, needs callable-mean quadrature)\n")
}

hdr("verdict")
tab    <- .results$tab
failed <- names(tab)[tab == "FAIL"]
cat(sprintf("  %d checks run, %d passed, %d failed.\n",
            length(tab), sum(tab == "PASS"), length(failed)))
if (length(failed)) {
  cat("  failing:\n"); for (f in failed) cat("    -", f, "\n")
  cat("  => do NOT write the extractor until these are understood.\n")
} else {
  cat("  => every alignment check passed: the in-scope extractor\n")
  cat("     (varying intercept on segment 1) is safe to write.\n")
}
cat("  Section C is descriptive, not a gate: it records what a\n")
cat("  callable-mean engine would have to handle for varying cps.\n")
