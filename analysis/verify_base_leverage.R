# =====================================================================
# Numerical scrutiny of the base-leverage decomposition (project: fibr_sbc)
#
# Claim to verify:
#   For the Fisher/observed-info metric G = [[G_BB,G_BF],[G_FB,G_FF]] of a
#   hierarchical posterior, base theta, fiber alpha, the per-observation
#   generalized (case-deletion) leverage splits ORTHOGONALLY:
#
#     g_i' G^{-1} g_i  =  g_i^F' G_FF^{-1} g_i^F   +   gtil_i' M^{-1} gtil_i
#                         \___ fiber (pooling) ___/    \___ base (horizontal) __/
#
#   gtil_i = g_i^B - G_BF G_FF^{-1} g_i^F = g_i^B + A' g_i^F ,  A = -G_FF^{-1} G_BF
#   M      = G_BB - G_BF G_FF^{-1} G_FB   (Schur complement = marginal base precision)
#
#   and the base term predicts each observation's LOO influence on sigma_u.
# =====================================================================

## ---- PART A: the identity is exact (random PD metric) ----------------
set.seed(1)
nB <- 3; nF <- 9
Z <- matrix(rnorm((nB+nF)^2), nB+nF); G <- crossprod(Z) + diag(nB+nF)  # random SPD
bi <- 1:nB; fi <- (nB+1):(nB+nF)
GBB<-G[bi,bi]; GBF<-G[bi,fi]; GFB<-G[fi,bi]; GFF<-G[fi,fi]
M <- GBB - GBF %*% solve(GFF) %*% GFB
g <- rnorm(nB+nF); gB<-g[bi]; gF<-g[fi]
gtil <- gB - GBF %*% solve(GFF) %*% gF
lhs <- as.numeric(t(g) %*% solve(G) %*% g)
rhs <- as.numeric(t(gF)%*%solve(GFF)%*%gF + t(gtil)%*%solve(M)%*%gtil)
cat(sprintf("PART A  identity: LHS=%.10f  RHS=%.10f  |diff|=%.2e  %s\n\n",
            lhs, rhs, abs(lhs-rhs), ifelse(abs(lhs-rhs)<1e-8,"EXACT","FAIL")))

## ---- PART B: Gaussian random-intercept LMM, validate base leverage ----
# y_ij ~ N(mu + a_j, sigma^2),  a_j ~ N(0, sigma_u^2).  Base=(mu, tau=log sigma_u),
# fiber = a.  sigma known.  A single observation's DIRECT score on tau is ZERO, so
# ALL of its influence on sigma_u must flow through the connection A -> cleanest
# possible test that the base leverage IS the connection channel.
set.seed(11)
sigma <- 1.0
J <- 12; n_j <- c(1,1,2,2,3,5,8,12,20,30,2,1); N <- sum(n_j)
mu_true <- 0; sigma_u_true <- 1.0
a_true <- rnorm(J, 0, sigma_u_true); a_true[3] <- 3.2; a_true[11] <- -3.0   # two outlier groups
grp <- rep(1:J, n_j)
y <- mu_true + a_true[grp] + rnorm(N, 0, sigma)

# mode given sigma_u (ridge / BLUP): solve linear system for (mu, a)
mode_given_su <- function(y, grp, sigma, sigma_u, keep = rep(TRUE, length(y))) {
  yy<-y[keep]; gg<-grp[keep]; nj<-tabulate(gg, J); ybar<-tapply(yy,factor(gg,1:J),sum)
  ybar[is.na(ybar)]<-0
  pu <- 1/sigma_u^2
  # normal equations for mu and a_j (a_j has prior precision pu)
  # d/da_j: (n_j/s2)(mu+a_j) - ybar_j/s2 + pu a_j = 0  -> a_j = (ybar_j/s2 - n_j mu/s2)/(n_j/s2+pu)
  # d/dmu:  sum_j n_j(mu+a_j) - sum y = 0
  s2<-sigma^2; w<-nj/s2; denom<-w+pu
  # a_j(mu) = (ybar_j/s2 - w_j mu)/denom_j ; plug into mu eqn
  A1 <- sum(nj) ; # coefficient bookkeeping
  # sum_j n_j (mu + a_j) = sum y  -> mu*sum n_j + sum n_j a_j = sum y
  # a_j = c_j - d_j mu ,  c_j=(ybar_j/s2)/denom_j, d_j=w_j/denom_j
  c_j <- (ybar/s2)/denom; d_j <- w/denom
  mu <- (sum(yy) - sum(nj*c_j)) / (sum(nj) - sum(nj*d_j))
  a <- c_j - d_j*mu
  list(mu=mu, a=as.numeric(a))
}
neg_log_post_tau <- function(tau, y, grp, sigma, keep=rep(TRUE,length(y))) {
  su<-exp(tau); m<-mode_given_su(y,grp,sigma,su,keep)
  yy<-y[keep]; gg<-grp[keep]; r<-yy-m$mu-m$a[gg]
  nj<-tabulate(gg,J)
  0.5*sum(r^2)/sigma^2 + 0.5*sum(m$a^2)/su^2 + sum((nj>0))*0 + J*tau + 0.5*tau^2  # tau prior N(0,1)
}
fit_full <- function(y,grp,sigma,keep=rep(TRUE,length(y))) {
  o<-optimize(neg_log_post_tau, c(-3,3), y=y, grp=grp, sigma=sigma, keep=keep)
  tau<-o$minimum; m<-mode_given_su(y,grp,sigma,exp(tau),keep); list(mu=m$mu, a=m$a, tau=tau)
}
fit <- fit_full(y,grp,sigma); mu<-fit$mu; a<-fit$a; tau<-fit$tau; su<-exp(tau)
cat(sprintf("PART B  fit: mu=%.3f  sigma_u=%.3f (tau=%.3f)\n", mu, su, tau))

# analytic G blocks at the mode (base=(mu,tau), fiber=a) -------------------
s2<-sigma^2; pu<-exp(-2*tau); nj<-tabulate(grp,J)
G_FF <- nj/s2 + pu                            # diagonal (vector)
Gmm <- N/s2 + 1/100                            # mu prior N(0,10^?) mild
Gtt <- 2*pu*sum(a^2) + 1                       # + tau prior precision 1
Gmt <- 0
G_BB <- matrix(c(Gmm,Gmt,Gmt,Gtt),2,2)
G_BF <- rbind(mu=nj/s2, tau=-2*a*pu)           # 2 x J
M <- G_BB - G_BF %*% (t(G_BF)/G_FF)            # 2x2 Schur (G_FF diagonal)
Minv <- solve(M)
pi_j <- pu/G_FF                                # fibr pooling factor per group

# per-observation leverages ------------------------------------------------
r <- y - mu - a[grp]                           # residuals
Lfib <- Lbase <- Ltot <- dtau_pred <- numeric(N)
for (i in 1:N) {
  j<-grp[i]; gF<-numeric(J); gF[j]<-r[i]/s2    # fiber score (only group j)
  gB<-c(mu=r[i]/s2, tau=0)                      # DIRECT base score: tau-part is ZERO
  GFFinv_gF <- gF/G_FF
  gtil <- gB - as.numeric(G_BF %*% GFFinv_gF)  # connection-corrected base score
  Lfib[i]  <- sum(gF^2/G_FF)
  Lbase[i] <- as.numeric(t(gtil)%*%Minv%*%gtil)
  Ltot[i]  <- Lfib[i]+Lbase[i]
  dtau_pred[i] <- -as.numeric(Minv %*% gtil)[2]  # one-step predicted LOO shift in tau
}

# ACTUAL leave-one-out shift in tau (refit without obs i) -------------------
dtau_act <- numeric(N)
for (i in 1:N) { keep<-rep(TRUE,N); keep[i]<-FALSE
  dtau_act[i] <- fit_full(y,grp,sigma,keep)$tau - tau }

cat(sprintf("PART B  cor(predicted dtau, actual dtau) = %.3f   [INFLUENCE -> sigma_u LOO shift]\n", cor(dtau_pred,dtau_act)))
cat(sprintf("PART B  cor(base influence, actual |dtau|) = %.3f\n", cor(Lbase, abs(dtau_act))))
# TWO distinct objects (the correction careful scrutiny forced):
#  - STRUCTURAL leverage h_i = info_i/G_FF,j : residual-free capacity; sum_j = 1-pi_j (= fibr).
#  - INFLUENCE  L_i = g_i'G^-1 g_i          : structural x realized standardized-residual^2; predicts THIS fold's k-hat.
h_struct   <- (1/s2)/G_FF[grp]                          # structural fiber leverage per obs
grp_struct <- as.numeric(tapply(h_struct, factor(grp,1:J), sum))
cat(sprintf("PART B  STRUCTURAL fiber leverage: max|sum_j h_i - (1-pi_j)| = %.2e  (identical to fibr)\n",
            max(abs(grp_struct - (1-pi_j)))))
grp_fiblev <- as.numeric(tapply(Lfib, factor(grp,1:J), sum))  # score-weighted (influence) version
cat(sprintf("PART B  (influence fiber term is residual-weighted, so cor with 1-pi_j is only %.2f -- expected)\n",
            cor(grp_fiblev, 1-pi_j)))

# ------------------------------- figure -----------------------------------
png("base_leverage_verify.png", width=1500, height=1050, res=140)
par(mfrow=c(2,3), mar=c(4.2,4.4,3,1))
# (1) predicted vs actual dtau  -- the validation
plot(dtau_pred, dtau_act, pch=19, col="#1f77b4", xlab="predicted LOO shift in tau (one-step)",
     ylab="actual LOO shift in tau (refit)", main=sprintf("Base leverage predicts sigma_u influence\ncor=%.3f", cor(dtau_pred,dtau_act)))
abline(0,1,col="grey50",lty=2)
# (2) base leverage vs group extremity, sized by n_j
plot(a[grp], Lbase, pch=19, col=ifelse(n_j[grp]<=2,"#d62728","#1f77b4"), cex=0.9,
     xlab=expression("group effect "*alpha[j]*" (obs's group)"), ylab="base leverage  L_base",
     main="Base leverage: extreme group x small n")
legend("top",bty="n",cex=0.85,pch=19,col=c("#d62728","#1f77b4"),legend=c("n_j <= 2 (near-singleton)","n_j > 2"))
# (3) STRUCTURAL fiber leverage == pooling factor (exact)
plot(1-pi_j, grp_struct, pch=19, col="#2ca02c",
     xlab=expression("1 - "*pi[j]*"  (fibr shrinkage)"), ylab="structural fiber leverage (summed)",
     main="Structural fiber leverage = 1 - pi_j (exact)"); abline(0,1,col="grey50",lty=2)
# (4) leverage decomposition stacked, ordered by total
o<-order(Ltot); barplot(rbind(Lfib[o],Lbase[o]), col=c("#2ca02c","#d62728"), border=NA,
     xlab="observations (ordered by total leverage)", ylab="generalized leverage",
     main="Total = fiber + base")
legend("topleft",bty="n",cex=0.85,fill=c("#2ca02c","#d62728"),legend=c("fiber (pooling)","base (horizontal)"))
# (5) which obs have high base leverage: the singletons / outliers
plot(1:N, Lbase, type="h", col="#d62728", xlab="observation index", ylab="base leverage",
     main="High base leverage = few-obs / outlier groups")
hi<-order(-Lbase)[1:4]; text(hi, Lbase[hi], labels=sprintf("g%d,n=%d",grp[hi],n_j[grp[hi]]), pos=3, cex=0.7)
# (6) base fraction of total leverage vs n_j
plot(n_j[grp], Lbase/Ltot, pch=19, col="#7f2704", log="x", xlab="group size n_j (log)",
     ylab="base share  L_base / L_total", main="Base share falls as groups grow")
dev.off()
cat("\nWrote figure: base_leverage_verify.png\n")
