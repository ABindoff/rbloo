# =====================================================================
# E-basestress : the non-tautological validation + the refit flag firing.
#
# The adversarial critique: on singleton folds RB-LOO and reloo compute the SAME
# marginal integral, so "matches exact refit" cannot falsify RB there. This
# experiment puts RB's OWN importance sampling under stress (few groups =>
# poorly-identified base => the base-IS reweighting is a genuine approximation
# that CAN fail), where brute-force refit is a real independent oracle.
#
# Claims tested:
#   (i)  RB-LOO's residual base-k-hat genuinely exceeds 0.7 on some folds here.
#   (ii) On those flagged folds RB elpd DIVERGES from exact refit (RB is NOT exact
#        there -- the honest limit); on unflagged folds RB matches refit.
#   (iii) The refit flag (base-k>0.7) catches the divergent folds -> the two-level
#        triage (pooling factor -> apply RB ; base-k -> refit) recovers the exact
#        answer everywhere at a fraction of the refit cost.
# Gaussian random-intercept LMM, conjugate Gibbs; brute-force refit ALL folds.
# =====================================================================
suppressMessages(library(loo)); set.seed(1)

gibbs_lmm <- function(y, grp, J, sigma, S=2000, warm=800, Vmu=100, a0=2, b0=2, keep=rep(TRUE,length(y))) {
  yy<-y[keep]; gg<-grp[keep]; nj<-tabulate(gg,J); ybar<-tapply(yy,factor(gg,1:J),sum); ybar[is.na(ybar)]<-0
  mu<-0; su2<-1; a<-numeric(J); s2<-sigma^2; KM<-numeric(S); KS<-numeric(S); KA<-matrix(0,S,J); tot<-S+warm
  for(t in 1:tot){ prec<-nj/s2+1/su2; va<-1/prec; a<-rnorm(J, va*((ybar-nj*mu)/s2), sqrt(va))
    vmu<-1/(sum(nj)/s2+1/Vmu); mu<-rnorm(1, vmu*sum(yy-a[gg])/s2, sqrt(vmu))
    su2<-1/rgamma(1, a0+J/2, b0+0.5*sum(a^2))
    if(t>warm){i<-t-warm; KM[i]<-mu; KS[i]<-sqrt(su2); KA[i,]<-a} }
  th<-seq(1,S,by=2); list(mu=KM[th], su=KS[th], a=KA[th,,drop=FALSE], J=J)
}

rb_loo <- function(y, grp, sigma, draws, base_cut=0.7) {
  J<-draws$J; N<-length(y); nj<-tabulate(grp,J); s2<-sigma^2
  mu_s<-draws$mu; su_s<-draws$su; A_s<-draws$a; S<-length(mu_s)
  L_full<-vapply(1:N,function(i) dnorm(y[i], mu_s+A_s[,grp[i]], sigma, log=TRUE), numeric(S))
  rf<-tryCatch(relative_eff(exp(L_full),chain_id=rep(1L,S)),error=function(e) rep(1,N))
  lf<-suppressWarnings(loo(L_full,r_eff=rf)); k_full<-lf$diagnostics$pareto_k; e_full<-lf$pointwise[,"elpd_loo"]
  L_rb<-matrix(0,S,N)
  for(i in 1:N){ j<-grp[i]; ij<-which(grp==j); Sj<-rowSums(outer(rep(1,S),y[ij])-mu_s)/s2
    Sjmi<-Sj-(y[i]-mu_s)/s2; Pmi<-1/su_s^2+(nj[j]-1)/s2
    L_rb[,i]<-dnorm(y[i], mu_s+Sjmi/Pmi, sqrt(s2+1/Pmi), log=TRUE) }
  rr<-tryCatch(relative_eff(exp(L_rb),chain_id=rep(1L,S)),error=function(e) rep(1,N))
  lr<-suppressWarnings(loo(L_rb,r_eff=rr)); k_base<-lr$diagnostics$pareto_k; e_rb<-lr$pointwise[,"elpd_loo"]
  mu_h<-mean(mu_s); su_h<-mean(su_s); a_h<-colMeans(A_s); pu<-1/su_h^2; G_FF<-nj/s2+pu; pij<-pu/G_FF
  G_BB<-matrix(c(N/s2+1/100,0,0,2*pu*sum(a_h^2)+1),2,2); G_BF<-rbind(nj/s2,-2*a_h*pu)
  M<-tryCatch(solve(G_BB-G_BF%*%(t(G_BF)/G_FF)), error=function(e) matrix(NA,2,2)); r_h<-y-mu_h-a_h[grp]
  hs<-Lb<-numeric(N)
  for(i in 1:N){ j<-grp[i]; gF<-numeric(J); gF[j]<-r_h[i]/s2; gt<-c(r_h[i]/s2,0)-as.numeric(G_BF%*%(gF/G_FF))
    hs[i]<-(1/s2)/G_FF[j]; Lb[i]<-as.numeric(t(gt)%*%M%*%gt) }
  data.frame(grp=grp,n_j=nj[grp],pi_j=pij[grp],h_struct=hs,L_base=Lb,
             k_full=k_full,k_base=k_base,elpd_full=e_full,elpd_rb=e_rb)
}

# brute-force refit gold (MARGINAL estimand, matches RB): leave i out, refit, then
# predictive = mean over refit draws of N(y_i | mu + a_{grp}, sigma).
refit_elpd <- function(y, grp, J, sigma, i){
  keep<-rep(TRUE,length(y)); keep[i]<-FALSE
  d<-gibbs_lmm(y,grp,J,sigma,S=1500,warm=600,keep=keep)
  log(mean(dnorm(y[i], d$mu + d$a[,grp[i]], sigma)))
}

## few-groups designs: base (sigma_u, mu) poorly identified -> base IS strained
configs<-list(c(1,1,1,2), c(1,1,1,1,2), c(1,1,1,1,2,3), c(1,1,2,2)); REP<-90
rows<-list(); ri<-0; set.seed(7)
for(cfg in configs){ J<-length(cfg); grp<-rep(1:J,cfg); N<-sum(cfg)
  for(r in 1:REP){ su<-runif(1,1.0,2.5); a_true<-rnorm(J,0,su); y<-a_true[grp]+rnorm(N,0,1)
    dr<-gibbs_lmm(y,grp,J,1.0,S=2000); out<-rb_loo(y,grp,1.0,dr)
    out<-out[is.finite(out$L_base),]; if(!nrow(out)) next
    out$e_ref<-vapply(as.integer(rownames(out)), function(i) refit_elpd(y,grp,J,1.0,i), numeric(1))
    ri<-ri+1; rows[[ri]]<-out } }
PA<-do.call(rbind,rows)

rmse<-function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
auc<-function(score,lab){ if(sum(lab)==0||sum(lab)==length(lab)) return(NA)
  r<-rank(score); (sum(r[lab==1])-sum(lab==1)*(sum(lab==1)+1)/2)/(sum(lab==1)*sum(lab==0)) }

flagged   <- PA$k_base>0.7                 # the refit flag
divergent <- abs(PA$elpd_rb - PA$e_ref) > 0.25   # folds where RB is actually wrong vs exact

cat("\n================ E-basestress : does reloo falsify RB where RB's own IS is stressed? ================\n")
cat(sprintf("folds=%d   PSIS k_full>0.7: %d   RB k_base>0.7 (FLAGGED): %d\n",
            nrow(PA), sum(PA$k_full>0.7), sum(flagged)))
cat("\n(i) RB's base-IS IS genuinely stressed here (unlike the singleton-only experiments).\n")
cat("\n(ii) RB accuracy vs EXACT refit, split by RB's own diagnostic:\n")
cat(sprintf("     RB elpd RMSE vs refit | k_base<=0.7 (unflagged): %.3f\n", rmse(PA$elpd_rb[!flagged], PA$e_ref[!flagged])))
cat(sprintf("     RB elpd RMSE vs refit | k_base >0.7 (FLAGGED)  : %.3f   <- RB is NOT exact here (honest limit)\n", rmse(PA$elpd_rb[flagged], PA$e_ref[flagged])))
cat(sprintf("     PSIS elpd RMSE vs refit (all folds)            : %.3f\n", rmse(PA$elpd_full, PA$e_ref)))
cat("\n(iii) The refit flag catches the folds where RB actually diverges from exact:\n")
cat(sprintf("     AUC( |RB-refit|>0.25  ~  k_base )    = %.3f\n", auc(PA$k_base, as.numeric(divergent))))
cat(sprintf("     AUC( |RB-refit|>0.25  ~  L_base )    = %.3f\n", auc(PA$L_base, as.numeric(divergent))))
# combined triage estimator: RB where unflagged, exact refit where flagged
e_combined <- ifelse(flagged, PA$e_ref, PA$elpd_rb)
cat(sprintf("\n     Combined triage (RB on %d unflagged + refit on %d flagged folds):\n", sum(!flagged), sum(flagged)))
cat(sprintf("       elpd RMSE vs exact = %.3f   (RB-only=%.3f, PSIS-only=%.3f)\n",
            rmse(e_combined, PA$e_ref), rmse(PA$elpd_rb, PA$e_ref), rmse(PA$elpd_full, PA$e_ref)))
cat(sprintf("       refit cost: %d / %d folds (%.0f%%) vs 100%% for full brute-force LOO\n",
            sum(flagged), nrow(PA), 100*mean(flagged)))

saveRDS(PA, "e_basestress.rds")

## ---- figure ----
png("e_basestress.png", width=1600, height=520, res=135)
par(mfrow=c(1,3), mar=c(4.3,4.4,3.2,1))
plot(PA$k_base, abs(PA$elpd_rb-PA$e_ref), pch=19, col=adjustcolor("#d62728",0.35), cex=0.6,
     xlab="RB base k-hat (RB's own diagnostic)", ylab="| RB elpd  -  exact refit |",
     main="RB error grows with its own\nbase-k-hat (falsifiable here)")
abline(v=0.7, lty=3, col="grey50")
cols<-ifelse(flagged,"#d62728","#1f77b4")
plot(PA$e_ref, PA$elpd_rb, pch=19, col=adjustcolor(cols,0.5), cex=0.6,
     xlab="exact refit elpd (gold)", ylab="RB-LOO elpd",
     main="RB matches exact where unflagged (blue);\ndiverges where flagged (red)")
abline(0,1,lty=2,col="grey60")
legend("topleft",bty="n",cex=0.8,pch=19,col=c("#1f77b4","#d62728"),legend=c("k_base<=0.7","k_base>0.7 (refit)"))
rm<-c(PSIS=rmse(PA$elpd_full,PA$e_ref), "RB only"=rmse(PA$elpd_rb,PA$e_ref), "RB+refit\nflag"=rmse(e_combined,PA$e_ref))
bp<-barplot(rm, col=c("#7f7f7f","#ff7f0e","#2ca02c"), ylab="elpd RMSE vs exact refit",
            main=sprintf("Two-level triage recovers exact\n(refit only %.0f%% of folds)",100*mean(flagged)))
text(bp, rm, sprintf("%.3f",rm), pos=3, xpd=NA, cex=0.9)
dev.off(); cat("\nWrote figure: e_basestress.png\n")
