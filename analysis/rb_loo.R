# =====================================================================
# rb_loo.R  --  P0 reference implementation + P1 gate-G1 experiment
# CV-reliability paper (project: fibr_sbc). Gaussian random-intercept LMM.
#   y_ij ~ N(mu + a_j, sigma^2),  a_j ~ N(0, sigma_u^2);  base=(mu, tau=log sigma_u).
# =====================================================================
suppressMessages(library(loo)); set.seed(1)

# ---- fit: conjugate Gibbs (sigma known) ----
gibbs_lmm <- function(y, grp, J, sigma, S=3000, warm=1000, Vmu=100, a0=2, b0=2,
                      keep=rep(TRUE,length(y))) {
  yy<-y[keep]; gg<-grp[keep]; nj<-tabulate(gg,J); ybar<-tapply(yy,factor(gg,1:J),sum); ybar[is.na(ybar)]<-0
  mu<-0; su2<-1; a<-numeric(J); s2<-sigma^2; KM<-numeric(S); KS<-numeric(S); KA<-matrix(0,S,J); tot<-S+warm
  for(t in 1:tot){
    prec<-nj/s2+1/su2; va<-1/prec; a<-rnorm(J, va*((ybar-nj*mu)/s2), sqrt(va))
    vmu<-1/(sum(nj)/s2+1/Vmu); mu<-rnorm(1, vmu*sum(yy-a[gg])/s2, sqrt(vmu))
    su2<-1/rgamma(1, a0+J/2, b0+0.5*sum(a^2))
    if(t>warm){i<-t-warm; KM[i]<-mu; KS[i]<-sqrt(su2); KA[i,]<-a}
  }
  th<-seq(1,S,by=2)                                   # light thinning -> ~independent for loo r_eff=1
  list(mu=KM[th], su=KS[th], a=KA[th,,drop=FALSE], J=J)
}

# =====================================================================
# rb_loo(): the three methods in one call.
#   draws = list(mu[S], su[S], a[S x J]) from the base+fiber posterior.
# Returns per-observation:
#   pi_j        Gelman-Pardoe pooling factor of the obs's group (fibr)
#   h_struct    structural fiber leverage (a-priori; sum over group = 1 - pi_j)
#   L_base      connection-generated base leverage (case-deletion influence on sigma_u)
#   k_full,elpd_full   full conditional PSIS-LOO (loo)          [the incumbent]
#   k_base,elpd_rb     Rao-Blackwellised LOO (fiber marginalised) [the cure]
#   refit_flag  TRUE where base IS is strained (k_base>0.7 OR L_base high)
# =====================================================================
rb_loo <- function(y, grp, sigma, draws, base_cut = 0.7) {
  J<-draws$J; N<-length(y); nj<-tabulate(grp,J); s2<-sigma^2
  mu_s<-draws$mu; su_s<-draws$su; A_s<-draws$a; S<-length(mu_s)

  ## full conditional PSIS-LOO
  L_full <- vapply(1:N, function(i) dnorm(y[i], mu_s + A_s[,grp[i]], sigma, log=TRUE), numeric(S))
  reff_f <- tryCatch(loo::relative_eff(exp(L_full), chain_id=rep(1L,S)), error=function(e) rep(1,N))
  lf <- suppressWarnings(loo::loo(L_full, r_eff=reff_f))
  k_full<-lf$diagnostics$pareto_k; elpd_full<-lf$pointwise[,"elpd_loo"]

  ## RB-LOO: marginalise a_j analytically given each base draw (rank-1 downdate)
  L_rb <- matrix(0,S,N)
  for(i in 1:N){ j<-grp[i]; ij<-which(grp==j)
    Sj   <- rowSums(outer(rep(1,S), y[ij]) - mu_s)/s2      # sum_{i' in j} (y-mu)/s2
    Sjmi <- Sj - (y[i]-mu_s)/s2                             # leave i out
    Pmi  <- 1/su_s^2 + (nj[j]-1)/s2                         # downdated group precision
    L_rb[,i] <- dnorm(y[i], mu_s + Sjmi/Pmi, sqrt(s2 + 1/Pmi), log=TRUE)
  }
  reff_r <- tryCatch(loo::relative_eff(exp(L_rb), chain_id=rep(1L,S)), error=function(e) rep(1,N))
  lr <- suppressWarnings(loo::loo(L_rb, r_eff=reff_r))
  k_base<-lr$diagnostics$pareto_k; elpd_rb<-lr$pointwise[,"elpd_loo"]

  ## analytic leverages at posterior mean (base=(mu,tau))
  mu_h<-mean(mu_s); su_h<-mean(su_s); a_h<-colMeans(A_s); pu<-1/su_h^2
  G_FF<-nj/s2+pu; pi_j_grp<-pu/G_FF
  G_BB<-matrix(c(N/s2+1/100,0,0,2*pu*sum(a_h^2)+1),2,2); G_BF<-rbind(nj/s2,-2*a_h*pu)
  M<-G_BB-G_BF%*%(t(G_BF)/G_FF); Minv<-solve(M); r_h<-y-mu_h-a_h[grp]
  h_struct<-L_base<-numeric(N)
  for(i in 1:N){ j<-grp[i]; gF<-numeric(J); gF[j]<-r_h[i]/s2
    gtil<-c(r_h[i]/s2,0)-as.numeric(G_BF%*%(gF/G_FF))
    h_struct[i]<-(1/s2)/G_FF[j]; L_base[i]<-as.numeric(t(gtil)%*%Minv%*%gtil) }

  data.frame(obs=1:N, grp=grp, n_j=nj[grp], pi_j=pi_j_grp[grp], h_struct=h_struct,
             L_base=L_base, k_full=k_full, k_base=k_base, elpd_full=elpd_full, elpd_rb=elpd_rb,
             refit_flag = (k_base>base_cut))
}

# =====================================================================
# P1: gate-G1 experiment. Replicate datasets; test H1, H2, H3 with MCSE.
# =====================================================================
sigma<-1.0; su_true<-1.3
n_j<-c(1,1,1,1,2,2,2,3,4,6,10,15); J<-length(n_j); N<-sum(n_j); grp<-rep(1:J,n_j)
REPS<-60; S_fit<-3000; nsub<-14

per_fold <- list(); acc <- list(); per_rep_rho_h<-numeric(REPS)
set.seed(100)
for(r in 1:REPS){
  a_true <- rnorm(J,0,su_true); y <- a_true[grp] + rnorm(N,0,sigma)
  dr <- gibbs_lmm(y,grp,J,sigma,S=S_fit)
  out <- rb_loo(y,grp,sigma,dr)
  # gold-standard refit on a leverage-stratified subset (for H2 accuracy)
  sub <- unique(c(order(-out$k_full)[1:6], order(-out$L_base)[1:4], order(out$k_full)[1:4]))[1:nsub]
  sub <- sub[!is.na(sub)]
  eref <- sapply(sub, function(i){ keep<-rep(TRUE,N); keep[i]<-FALSE
    d<-gibbs_lmm(y,grp,J,sigma,S=1500,warm=500,keep=keep)
    log(mean(dnorm(y[i], d$mu+d$a[,grp[i]], sigma))) })
  out$rep<-r; per_fold[[r]]<-out
  acc[[r]]<-data.frame(rep=r, obs=sub, k_full=out$k_full[sub],
                       e_full=out$elpd_full[sub], e_rb=out$elpd_rb[sub], e_ref=eref)
  per_rep_rho_h[r]<-suppressWarnings(cor(out$h_struct,out$k_full,method="spearman"))
  if(r%%15==0) cat(sprintf("  ...replicate %d/%d\n", r, REPS))
}
PF<-do.call(rbind,per_fold); AC<-do.call(rbind,acc)

## ---- H1: pooling/structural leverage predicts k_full (a-priori) ----
rho_pool <- suppressWarnings(cor(PF$h_struct, PF$k_full, method="spearman"))
mcse_rho <- sd(per_rep_rho_h)/sqrt(REPS)
lab <- as.numeric(PF$k_full>0.7)
auc <- function(score,lab){ if(sum(lab)==0||sum(lab)==length(lab)) return(NA)
  r<-rank(score); (sum(r[lab==1])-sum(lab==1)*(sum(lab==1)+1)/2)/(sum(lab==1)*sum(lab==0)) }
auc_pool <- auc(PF$h_struct, lab)

## ---- H2: RB-LOO vs PSIS-LOO accuracy against refit ----
rmse<-function(a,b) sqrt(mean((a-b)^2))
rmse_rb<-rmse(AC$e_rb,AC$e_ref); rmse_ps<-rmse(AC$e_full,AC$e_ref)
rb_closer <- mean(abs(AC$e_rb-AC$e_ref) < abs(AC$e_full-AC$e_ref))
# focus on high-k folds
hi<-AC$k_full>0.7
rmse_rb_hi<-if(any(hi)) rmse(AC$e_rb[hi],AC$e_ref[hi]) else NA
rmse_ps_hi<-if(any(hi)) rmse(AC$e_full[hi],AC$e_ref[hi]) else NA

## ---- H3: RB cures fiber failures ----
badf<-PF$k_full>0.7; cured<-mean(PF$k_base[badf]<0.7)
mean_kfull<-mean(PF$k_full); mean_kbase<-mean(PF$k_base)

cat("\n================ GATE G1 RESULTS ================\n")
cat(sprintf("replicates=%d  folds=%d  (k_full>0.7: %d folds; refit subset n=%d)\n",
            REPS, nrow(PF), sum(badf), nrow(AC)))
cat(sprintf("H1  Spearman(pooling/structural leverage, k_full) = %.3f  (MCSE %.3f)   [thr >=0.60 %s]\n",
            rho_pool, mcse_rho, ifelse(rho_pool>=0.6,"PASS","FAIL")))
cat(sprintf("H1  AUC(k_full>0.7 ~ structural leverage)         = %.3f            [thr >=0.80 %s]\n",
            auc_pool, ifelse(!is.na(auc_pool)&&auc_pool>=0.8,"PASS","FAIL")))
cat(sprintf("H2  elpd RMSE vs refit:  RB-LOO=%.4f  PSIS-LOO=%.4f  ratio=%.2f  [thr <=0.50 %s]\n",
            rmse_rb, rmse_ps, rmse_rb/rmse_ps, ifelse(rmse_rb/rmse_ps<=0.5,"PASS","FAIL")))
cat(sprintf("H2  RB closer to refit than PSIS: %.0f%% of subset folds\n", 100*rb_closer))
cat(sprintf("H2  high-k folds only: RB RMSE=%.4f  PSIS RMSE=%.4f\n", rmse_rb_hi, rmse_ps_hi))
cat(sprintf("H3  cured (k_base<0.7 | k_full>0.7) = %.0f%%   (mean k: %.2f -> %.2f)  [thr >=80%% %s]\n",
            100*cured, mean_kfull, mean_kbase, ifelse(cured>=0.8,"PASS","FAIL")))

# ---- figure ----
png("rb_loo_gateG1.png", width=1550, height=520, res=135)
par(mfrow=c(1,3), mar=c(4.2,4.3,3,1))
cc<-densCols <- adjustcolor(ifelse(PF$k_base>0.7,"#d62728","#1f77b4"),0.35)
plot(PF$h_struct, PF$k_full, pch=19, col=adjustcolor("#333333",0.25), cex=0.6,
     xlab="pooling / structural leverage (a-priori)", ylab="PSIS-LOO k-hat (full)",
     main=sprintf("H1: predicts k-hat\nrho=%.2f, AUC=%.2f", rho_pool, auc_pool)); abline(h=0.7,lty=3,col="red")
plot(PF$k_full, PF$k_base, pch=19, col=adjustcolor("#2ca02c",0.25), cex=0.6,
     xlab="k-hat full (conditional)", ylab="k-hat base (RB-LOO)",
     main=sprintf("H3: RB cures failures\n%.0f%% of k>0.7 folds cured", 100*cured))
abline(0,1,lty=2,col="grey60"); abline(h=0.7,v=0.7,lty=3,col="grey70")
plot(AC$e_ref, AC$e_rb, pch=19, col="#1f77b4", cex=0.7, xlab="refit elpd (gold)", ylab="LOO elpd",
     main=sprintf("H2: accuracy\nRMSE RB=%.3f  PSIS=%.3f", rmse_rb, rmse_ps))
points(AC$e_ref, AC$e_full, pch=1, col="#d62728", cex=0.7); abline(0,1,lty=2,col="grey60")
legend("topleft",bty="n",cex=0.85,pch=c(19,1),col=c("#1f77b4","#d62728"),legend=c("RB-LOO","PSIS-LOO"))
dev.off(); cat("\nWrote figure: rb_loo_gateG1.png\n")

saveRDS(list(PF=PF,AC=AC,per_rep_rho_h=per_rep_rho_h), "rb_loo_gateG1.rds")
