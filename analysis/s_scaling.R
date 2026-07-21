# =====================================================================
# E-Sscaling : accuracy vs number of posterior draws S.
# The "efficient" claim, made concrete: on the failing folds, RB-LOO's weights
# have finite variance so its elpd converges fast in S, while PSIS-LOO's
# heavy-tailed weights leave it noisy and biased -- needing far more draws and
# never fully catching up. Gaussian random-intercept LMM (conjugate Gibbs); the
# per-fold estimator is subsampled from one long chain so the only thing varying
# is S. Gold = high-precision brute-force refit.
# =====================================================================
suppressMessages(library(loo)); set.seed(1)

gibbs_lmm <- function(y, grp, J, sigma, S=8000, warm=1500, Vmu=100, a0=2, b0=2, keep=rep(TRUE,length(y))) {
  yy<-y[keep]; gg<-grp[keep]; nj<-tabulate(gg,J); ybar<-tapply(yy,factor(gg,1:J),sum); ybar[is.na(ybar)]<-0
  mu<-0; su2<-1; a<-numeric(J); s2<-sigma^2; KM<-numeric(S); KS<-numeric(S); KA<-matrix(0,S,J); tot<-S+warm
  for(t in 1:tot){ prec<-nj/s2+1/su2; va<-1/prec; a<-rnorm(J, va*((ybar-nj*mu)/s2), sqrt(va))
    vmu<-1/(sum(nj)/s2+1/Vmu); mu<-rnorm(1, vmu*sum(yy-a[gg])/s2, sqrt(vmu))
    su2<-1/rgamma(1, a0+J/2, b0+0.5*sum(a^2))
    if(t>warm){i<-t-warm; KM[i]<-mu; KS[i]<-sqrt(su2); KA[i,]<-a} }
  list(mu=KM, su=KS, a=KA, J=J)
}
# per-fold PSIS-LOO and RB-LOO elpd from a draws subset
elpds <- function(y, grp, sigma, dr) {
  J<-dr$J; N<-length(y); nj<-tabulate(grp,J); s2<-sigma^2
  mu_s<-dr$mu; su_s<-dr$su; A_s<-dr$a; S<-length(mu_s)
  L_full<-vapply(1:N,function(i) dnorm(y[i], mu_s+A_s[,grp[i]], sigma, log=TRUE), numeric(S))
  rf<-tryCatch(relative_eff(exp(L_full),chain_id=rep(1L,S)),error=function(e) rep(1,N))
  lf<-suppressWarnings(loo(L_full,r_eff=rf))
  L_rb<-matrix(0,S,N)
  for(i in 1:N){ j<-grp[i]; ij<-which(grp==j); Sj<-rowSums(outer(rep(1,S),y[ij])-mu_s)/s2
    Sjmi<-Sj-(y[i]-mu_s)/s2; Pmi<-1/su_s^2+(nj[j]-1)/s2
    L_rb[,i]<-dnorm(y[i], mu_s+Sjmi/Pmi, sqrt(s2+1/Pmi), log=TRUE) }
  rr<-tryCatch(relative_eff(exp(L_rb),chain_id=rep(1L,S)),error=function(e) rep(1,N))
  lr<-suppressWarnings(loo(L_rb,r_eff=rr))
  list(k=lf$diagnostics$pareto_k, e_full=lf$pointwise[,"elpd_loo"], e_rb=lr$pointwise[,"elpd_loo"])
}
refit_gold <- function(y, grp, J, sigma, i){
  keep<-rep(TRUE,length(y)); keep[i]<-FALSE
  d<-gibbs_lmm(y,grp,J,sigma,S=6000,warm=1500,keep=keep)
  log(mean(dnorm(y[i], d$mu + d$a[,grp[i]], sigma)))
}

n_j<-c(1,1,1,1,2,2,3,5,8); J<-length(n_j); N<-sum(n_j); grp<-rep(1:J,n_j)
sigma<-1.0; su_true<-1.4; S_max<-8000; S_grid<-c(250,500,1000,2000,4000,8000)
REPS<-12
err_ps<-err_rb<-matrix(0,length(S_grid),0)   # cols = folds pooled across reps (high-k)
set.seed(50)
for(r in 1:REPS){
  a_true<-rnorm(J,0,su_true); y<-a_true[grp]+rnorm(N,0,sigma)
  dr<-gibbs_lmm(y,grp,J,sigma,S=S_max)
  ref<-elpds(y,grp,sigma,dr)                      # full-S reference to pick high-k folds
  hi<-which(ref$k>0.7); if(!length(hi)) next
  gold<-vapply(hi, function(i) refit_gold(y,grp,J,sigma,i), numeric(1))
  Ep<-Er<-matrix(NA,length(S_grid),length(hi))
  for(si in seq_along(S_grid)){ S<-S_grid[si]
    drs<-list(mu=dr$mu[1:S], su=dr$su[1:S], a=dr$a[1:S,,drop=FALSE], J=J)
    e<-elpds(y,grp,sigma,drs)
    Ep[si,]<-e$e_full[hi]-gold; Er[si,]<-e$e_rb[hi]-gold }
  err_ps<-cbind(err_ps,Ep); err_rb<-cbind(err_rb,Er)
  cat(sprintf("  rep %2d/%d : %d high-k folds (cumulative %d)\n", r,REPS,length(hi),ncol(err_ps)))
}
rmse_ps<-sqrt(rowMeans(err_ps^2)); rmse_rb<-sqrt(rowMeans(err_rb^2))
cat("\n================ E-Sscaling : elpd RMSE vs refit on high-k folds ================\n")
cat(sprintf("  pooled high-k folds: %d\n", ncol(err_ps)))
for(si in seq_along(S_grid)) cat(sprintf("  S=%5d :  PSIS-LOO RMSE=%.3f   RB-LOO RMSE=%.3f\n", S_grid[si], rmse_ps[si], rmse_rb[si]))
saveRDS(list(S_grid=S_grid, rmse_ps=rmse_ps, rmse_rb=rmse_rb, err_ps=err_ps, err_rb=err_rb), "s_scaling.rds")

png("s_scaling.png", width=900, height=560, res=135)
par(mar=c(4.6,4.6,3,1))
yl<-range(c(rmse_ps,rmse_rb)); plot(S_grid, rmse_ps, type="b", pch=1, col="#d62728", lwd=2, log="x",
     ylim=yl, xlab="posterior draws S (log scale)", ylab="elpd RMSE vs exact refit (high-k folds)",
     main="RB-LOO is accurate at small S;\nPSIS-LOO needs more draws and plateaus above it")
lines(S_grid, rmse_rb, type="b", pch=19, col="#1f77b4", lwd=2)
legend("topright", bty="n", pch=c(1,19), col=c("#d62728","#1f77b4"), lwd=2, legend=c("PSIS-LOO","RB-LOO"))
dev.off(); cat("\nWrote figure: s_scaling.png\n")
