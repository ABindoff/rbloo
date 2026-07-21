# =====================================================================
# Rao-Blackwellised PIT for Simulation-Based Calibration  --  Phase 1 prototype
# Project: fibr_sbc
#
# Model (conjugate normal-normal hierarchy, sigma_y known):
#     y_ij ~ N(alpha_j, sigma_y^2),   i = 1..n_j
#     alpha_j ~ N(mu, sigma^2)                 (fiber / group effects)
#     mu ~ N(0, tau^2)                         (base)
#     sigma^2 ~ InvGamma(a0, b0)               (base)
#
# We compare two SBC rank statistics computed from the SAME correct-or-buggy
# Gibbs chain:
#   (1) counting rank   u_cnt = mean_l 1{alpha_j^(l) < alpha_j~}      (standard SBC)
#   (2) Rao-Blackwellised PIT
#         u_rb = mean_l  Phi( (alpha_j~ - m_j^(l)) / s_j^(l) )
#       where, GIVEN the base draw (mu^(l), sigma^(l)), the fiber conditional is
#       exactly Gaussian with
#         prec_j^(l) = 1/sigma^(l)^2 + n_j/sigma_y^2     = G_FF,j   (fibr's total precision)
#         s_j^(l)    = 1/sqrt(prec_j^(l))
#         pi_j^(l)   = (1/sigma^(l)^2) / prec_j^(l)       = fibr prior fraction
#         m_j^(l)    = pi_j^(l)*mu^(l) + (1 - pi_j^(l))*ybar_j
#
#   The pi_j / prec_j pair is EXACTLY what fibr::prior_fraction() returns
#   (prior precision, lik_information = n_j/sigma_y^2). In production the block
#   below is replaced by a call to fibr. Here we inline the identical formula so
#   the prototype is dependency-free.
# =====================================================================

set.seed(1)

# ---- fibr prior-fraction formula (faithful copy of fibr::prior_fraction.default) ----
# returns list(pi, prec, s) given per-draw prior precision and fixed lik information
fibr_moments <- function(prior_prec, lik_info) {
  prec <- prior_prec + lik_info
  list(pi = prior_prec / prec, prec = prec, s = 1 / sqrt(prec))
}

# ---------------------------- model settings -------------------------------
J        <- 6
n_j      <- c(2, 5, 10, 20, 40, 80)     # spread of group sizes -> spread of prior fractions
sigma_y  <- 1.0
tau      <- 5.0                          # prior sd on mu
a0       <- 3.0; b0 <- 2.0               # InvGamma(sigma^2)

# ------------------------- prior-predictive draw ---------------------------
sim_prior <- function() {
  mu     <- rnorm(1, 0, tau)
  sigma2 <- 1 / rgamma(1, shape = a0, rate = b0)
  alpha  <- rnorm(J, mu, sqrt(sigma2))
  ybar   <- rnorm(J, alpha, sigma_y / sqrt(n_j))   # sufficient stat: group means
  sse    <- sum((n_j - 1))                          # not needed (sigma_y known)
  list(mu = mu, sigma2 = sigma2, alpha = alpha, ybar = ybar)
}

# ------------------------------- Gibbs -------------------------------------
# mu_infl / alpha_infl < 1  ==  posterior conditional made too NARROW (a bug).
#   mu_infl    : bug in the BASE (mu) update
#   alpha_infl : bug in the FIBER (alpha) update  -- the coordinate we RB away
#   s2_bias    : bug in the variance-component (sigma^2) update -- a systematic
#                scale error, propagates into every group's fiber conditional.
gibbs <- function(data, L = 100, warmup = 50, alpha_infl = 1, s2_bias = 1) {
  ybar <- data$ybar
  mu <- 0; sigma2 <- 1
  keep_mu <- numeric(L); keep_s2 <- numeric(L)
  keep_al <- matrix(0, L, J)
  tot <- warmup + L
  for (t in 1:tot) {
    # alpha | .
    prec_a <- 1 / sigma2 + n_j / sigma_y^2
    v_a    <- 1 / prec_a
    m_a    <- v_a * (mu / sigma2 + n_j * ybar / sigma_y^2)
    alpha  <- rnorm(J, m_a, sqrt(v_a) * alpha_infl)
    # mu | .
    v_mu   <- 1 / (J / sigma2 + 1 / tau^2)
    m_mu   <- v_mu * (sum(alpha) / sigma2)
    mu     <- rnorm(1, m_mu, sqrt(v_mu))
    # sigma^2 | .   (s2_bias != 1 injects a systematic scale error)
    sigma2 <- s2_bias / rgamma(1, shape = a0 + J / 2, rate = b0 + 0.5 * sum((alpha - mu)^2))
    if (t > warmup) { i <- t - warmup; keep_mu[i] <- mu; keep_s2[i] <- sigma2; keep_al[i, ] <- alpha }
  }
  list(mu = keep_mu, sigma2 = keep_s2, alpha = keep_al)
}

# --------------------------- the two PITs ----------------------------------
pit_counting <- function(chain, alpha_true) {
  vapply(1:J, function(j) mean(chain$alpha[, j] < alpha_true[j]), numeric(1))
}
pit_rb <- function(chain, data) {   # Rao-Blackwellised over the base draws
  prior_prec <- 1 / chain$sigma2               # per-draw
  vapply(1:J, function(j) {
    mo <- fibr_moments(prior_prec, n_j[j] / sigma_y^2)   # <- fibr call in production
    m  <- mo$pi * chain$mu + (1 - mo$pi) * data$ybar[j]
    mean(pnorm((data$alpha[j] - m) / mo$s))
  }, numeric(1))
}

# ---------------- simultaneous ECDF band (KS-type, MC null) ----------------
ks_crit <- function(N, level = 0.05, M = 4000) {   # continuous-uniform null (for RB PIT)
  D <- numeric(M)
  for (m in 1:M) {
    u <- sort(runif(N)); i <- 1:N
    D[m] <- max(pmax(i / N - u, u - (i - 1) / N))
  }
  as.numeric(quantile(D, 1 - level))
}
# Discrete-uniform null for the COUNTING rank at L draws: under calibration the
# rank is uniform on {0,..,L}, so the PIT lives on {0,1/L,..,1}. Using this null
# (not the continuous one) is what makes small-L counting ranks a fair comparison.
ks_crit_discrete <- function(N, Ldraws, level = 0.05, M = 4000) {
  D <- numeric(M)
  for (m in 1:M) {
    u <- sample(0:Ldraws, N, replace = TRUE) / Ldraws
    D[m] <- ks_stat(u)
  }
  as.numeric(quantile(D, 1 - level))
}
ks_stat <- function(u) {
  u <- sort(u); N <- length(u); i <- 1:N
  max(pmax(i / N - u, u - (i - 1) / N))
}

# ============================ EXPERIMENT 1 =================================
# Correctness: correct sampler -> both PITs uniform (inside simultaneous band)
run_sbc <- function(S, L, s2_bias = 1, alpha_infl = 1) {
  cnt <- matrix(0, S, J); rb <- matrix(0, S, J)
  for (s in 1:S) {
    d  <- sim_prior()
    ch <- gibbs(d, L = L, s2_bias = s2_bias, alpha_infl = alpha_infl)
    cnt[s, ] <- pit_counting(ch, d$alpha)
    rb[s, ]  <- pit_rb(ch, d)
  }
  list(cnt = as.vector(cnt), rb = as.vector(rb))
}

cat("Experiment 1: correctness (correct sampler)\n")
S1 <- 400; L1 <- 150
e1 <- run_sbc(S1, L1)
crit1 <- ks_crit(length(e1$cnt))
cat(sprintf("  KS stat  counting = %.4f | RB = %.4f | 0.05 crit = %.4f (N=%d)\n",
            ks_stat(e1$cnt), ks_stat(e1$rb), crit1, length(e1$cnt)))

# ============================ EXPERIMENT 2 =================================
# Efficiency: fix ONE dataset, rerun the sampler many times, measure the
# within-simulation SD of each PIT estimate. RB integrates alpha out analytically
# so only base-draw noise remains -> much smaller SD -> effective-L multiplier.
cat("Experiment 2: within-simulation variance reduction\n")
d_fixed <- sim_prior()
R_rep <- 300; Lfix <- 100
cnt_rep <- matrix(0, R_rep, J); rb_rep <- matrix(0, R_rep, J)
for (r in 1:R_rep) {
  ch <- gibbs(d_fixed, L = Lfix)
  cnt_rep[r, ] <- pit_counting(ch, d_fixed$alpha)
  rb_rep[r, ]  <- pit_rb(ch, d_fixed)
}
sd_cnt <- apply(cnt_rep, 2, sd); sd_rb <- apply(rb_rep, 2, sd)
eff_L  <- (sd_cnt / sd_rb)^2      # L would need to grow by this factor to match RB
pf_fixed <- fibr_moments(1 / mean(gibbs(d_fixed, L = 200)$sigma2), n_j / sigma_y^2)$pi
for (j in 1:J)
  cat(sprintf("  group %d (n=%2d, pi~%.2f): SD_cnt=%.4f SD_rb=%.4f  effective-L x%.1f\n",
              j, n_j[j], pf_fixed[j], sd_cnt[j], sd_rb[j], eff_L[j]))

# ============================ EXPERIMENT 3 =================================
# Power vs L (posterior draws per fit) under a BASE bug (variance component
# sigma^2 biased low). This is the core speed lever: the counting rank needs many
# draws L to resolve the miscalibration; the RB PIT is near-exact for any L, so it
# holds full power as L -> a handful. Fixed S; sweep L.
cat("Experiment 3: power vs L under a base (sigma^2) bug (fixed S); fair per-method nulls\n")
Lvec <- c(3, 5, 10, 25, 60, 150); B <- 40; S3 <- 400; s2_bug <- 0.70
crit_rb3 <- ks_crit(S3 * J, M = 3000)               # continuous null for RB
pow_cnt <- pow_rb <- numeric(length(Lvec))
size_cnt <- numeric(length(Lvec))                    # type-I check for counting
for (k in seq_along(Lvec)) {
  L <- Lvec[k]
  crit_cnt <- ks_crit_discrete(S3 * J, L, M = 3000)  # discrete null for counting at this L
  rej_c <- rej_r <- fp_c <- 0
  for (b in 1:B) {
    e  <- run_sbc(S3, L, s2_bias = s2_bug)           # buggy sampler
    e0 <- run_sbc(S3, L)                             # correct sampler (size check)
    rej_c <- rej_c + (ks_stat(e$cnt) > crit_cnt)
    rej_r <- rej_r + (ks_stat(e$rb)  > crit_rb3)
    fp_c  <- fp_c  + (ks_stat(e0$cnt) > crit_cnt)
  }
  pow_cnt[k] <- rej_c / B; pow_rb[k] <- rej_r / B; size_cnt[k] <- fp_c / B
  cat(sprintf("  L=%3d : power counting=%.2f  RB=%.2f  (counting size=%.2f)\n",
              L, pow_cnt[k], pow_rb[k], size_cnt[k]))
}

# ============================ EXPERIMENT 4 =================================
# The honest caveat: a FIBER bug (alpha update too narrow). Counting ranks see
# it; the RB PIT -- which replaces the sampler's alpha draws with the correct
# analytic conditional -- is BLIND to it. RB relocates what is under test.
cat("Experiment 4: fiber (alpha) bug -- counting detects, RB blind\n")
S4 <- 400; L4 <- 150; a_bug <- 0.80
e4 <- run_sbc(S4, L4, alpha_infl = a_bug)
crit4 <- ks_crit(length(e4$cnt))
cat(sprintf("  KS stat  counting = %.4f | RB = %.4f | 0.05 crit = %.4f\n",
            ks_stat(e4$cnt), ks_stat(e4$rb), crit4))

# =============================== FIGURE ====================================
ecdf_diff <- function(u) { u <- sort(u); list(x = u, y = (1:length(u))/length(u) - u) }
panel <- function(u, crit, main, col) {
  ed <- ecdf_diff(u)
  plot(ed$x, ed$y, type = "s", col = col, lwd = 1.6, ylim = c(-1,1)*max(crit*1.6, max(abs(ed$y))),
       xlab = "PIT value", ylab = expression(hat(F)(u) - u), main = main, cex.main = 0.98)
  abline(h = 0, col = "grey70")
  polygon(c(0,1,1,0), c(crit,crit,-crit,-crit), col = adjustcolor("grey50",0.14), border = NA)
  inside <- ks_stat(u) <= crit
  legend("topright", bty = "n", cex = 0.85,
         legend = sprintf("D=%.3f (%s)", ks_stat(u), ifelse(inside,"in band","OUT")))
}

png("rb_pit_sbc.png", width = 1500, height = 1450, res = 140)
par(mfrow = c(3,2), mar = c(4,4.3,3,1))
col_c <- "#1f77b4"; col_r <- "#d62728"
panel(e1$cnt, crit1, "E1  Correct sampler - counting rank", col_c)
panel(e1$rb,  crit1, "E1  Correct sampler - Rao-Blackwellised PIT", col_r)

# power curve vs L
plot(Lvec, pow_cnt, type = "b", pch = 19, col = col_c, ylim = c(0,1), lwd = 2, log = "x",
     xlab = "posterior draws per fit  L  (log)", ylab = "detection power (alpha=0.05)",
     main = expression("E3  Power vs L under base ("*sigma^2*") bug, S=400"))
lines(Lvec, pow_rb, type = "b", pch = 17, col = col_r, lwd = 2)
abline(h = 0.05, lty = 3, col = "grey60")
legend("right", bty = "n", legend = c("counting rank","Rao-Blackwellised"),
       col = c(col_c,col_r), pch = c(19,17), lwd = 2, cex = 0.9)

# efficiency bars
bp <- barplot(rbind(sd_cnt, sd_rb), beside = TRUE, col = c(col_c,col_r),
              names.arg = paste0("n=", n_j), ylab = "within-sim SD of PIT",
              main = "E2  Variance reduction (fixed dataset)", cex.names = 0.8,
              ylim = c(0, 0.062))
legend("top", bty = "n", legend = c("counting","RB (analytic fiber)"),
       fill = c(col_c,col_r), cex = 0.85, horiz = TRUE)
text(colMeans(bp), pmax(sd_cnt, sd_rb) + 0.003, sprintf("x%.0f", eff_L), cex = 0.8, col = "grey20")

panel(e4$cnt, crit4, "E4  Fiber (alpha) bug - counting rank DETECTS", col_c)
panel(e4$rb,  crit4, "E4  Fiber (alpha) bug - RB PIT is BLIND", col_r)
dev.off()
cat("\nWrote figure: rb_pit_sbc.png\n")
