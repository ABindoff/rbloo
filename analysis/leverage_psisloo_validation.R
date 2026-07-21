# =====================================================================
# Validation Model 1: leverage decomposition vs PSIS-LOO k-hat
# Gaussian random-intercept LMM, conjugate Gibbs (no Stan). (project: fibr_sbc)
# Two regimes: A = sigma_u well identified (15 groups); B = barely (6 groups).
# See EXPERIMENT_leverage_validation.md for the pre-registered hypotheses.
# =====================================================================
suppressMessages(library(loo)); set.seed(20)
sigma <- 1.0

gibbs <- function(y, grp, J, sigma, S=3000, warm=1000, Vmu=100, a0=2, b0=2, keep=rep(TRUE,length(y))) {
  yy<-y[keep]; gg<-grp[keep]; nj<-tabulate(gg,J); ybar<-tapply(yy,factor(gg,1:J),sum); ybar[is.na(ybar)]<-0
  mu<-0; su2<-1; a<-numeric(J); s2<-sigma^2; KM<-numeric(S); KS<-numeric(S); KA<-matrix(0,S,J); tot<-S+warm
  for(t in 1:tot){
    prec<-nj/s2+1/su2; va<-1/prec; a<-rnorm(J, va*((ybar-nj*mu)/s2), sqrt(va))
    vmu<-1/(sum(nj)/s2+1/Vmu); mu<-rnorm(1, vmu*sum(yy-a[gg])/s2, sqrt(vmu))
    su2<-1/rgamma(1, a0+J/2, b0+0.5*sum(a^2))
    if(t>warm){i<-t-warm; KM[i]<-mu; KS[i]<-sqrt(su2); KA[i,]<-a}
  }
  list(mu=KM, su=KS, a=KA)
}
sp<-function(a,b) suppressWarnings(cor(a,b,method="spearman"))

run_scenario <- function(n_j, a_true, label) {
  J<-length(n_j); N<-sum(n_j); grp<-rep(1:J,n_j); s2<-sigma^2
  y <- a_true[grp] + rnorm(N,0,sigma)
  dr<-gibbs(y,grp,J,sigma); idx<-seq(1,length(dr$mu),by=2)
  mu_s<-dr$mu[idx]; su_s<-dr$su[idx]; A_s<-dr$a[idx,,drop=FALSE]; S<-length(mu_s); nj<-tabulate(grp,J)
  # full conditional PSIS-LOO
  L_full<-sapply(1:N,function(i) dnorm(y[i], mu_s+A_s[,grp[i]], sigma, log=TRUE))
  lf<-suppressWarnings(loo(L_full,r_eff=rep(1,N))); k_full<-lf$diagnostics$pareto_k; e_full<-lf$pointwise[,"elpd_loo"]
  # RB-LOO: fiber marginalised (analytic downdate) given base draws
  L_rb<-matrix(0,S,N)
  for(i in 1:N){ j<-grp[i]; ij<-which(grp==j)
    Sj<-rowSums(outer(rep(1,S),y[ij])-mu_s)/s2; Sj_mi<-Sj-(y[i]-mu_s)/s2
    P_mi<-1/su_s^2+(nj[j]-1)/s2; L_rb[,i]<-dnorm(y[i], mu_s+Sj_mi/P_mi, sqrt(s2+1/P_mi), log=TRUE) }
  lr<-suppressWarnings(loo(L_rb,r_eff=rep(1,N))); k_base<-lr$diagnostics$pareto_k; e_rb<-lr$pointwise[,"elpd_loo"]
  # analytic leverages at posterior mean
  mu_h<-mean(mu_s); su_h<-mean(su_s); a_h<-colMeans(A_s); pu<-1/su_h^2; G_FF<-nj/s2+pu; pi_j<-pu/G_FF
  G_BB<-matrix(c(N/s2+1/100,0,0,2*pu*sum(a_h^2)+1),2,2); G_BF<-rbind(nj/s2,-2*a_h*pu)
  M<-G_BB-G_BF%*%(t(G_BF)/G_FF); Minv<-solve(M); r_h<-y-mu_h-a_h[grp]
  Lfib<-Lbase<-Ltot<-hstruct<-numeric(N)
  for(i in 1:N){ j<-grp[i]; gF<-numeric(J); gF[j]<-r_h[i]/s2; gtil<-c(r_h[i]/s2,0)-as.numeric(G_BF%*%(gF/G_FF))
    Lfib[i]<-sum(gF^2/G_FF); Lbase[i]<-as.numeric(t(gtil)%*%Minv%*%gtil); Ltot[i]<-Lfib[i]+Lbase[i]; hstruct[i]<-(1/s2)/G_FF[j] }
  # refit gold standard (stratified subset)
  sub<-unique(c(order(-k_full)[1:6], order(-Lbase)[1:4], order(k_full)[1:4], if(N>14) sample(N,3) else integer(0)))
  e_ref<-sapply(sub,function(i){ keep<-rep(TRUE,N); keep[i]<-FALSE; d<-gibbs(y,grp,J,sigma,S=1500,warm=500,keep=keep)
    log(mean(dnorm(y[i], d$mu+d$a[,grp[i]], sigma))) })
  list(label=label,N=N,J=J,n_j=n_j,grp=grp,k_full=k_full,k_base=k_base,e_full=e_full,e_rb=e_rb,
       Ltot=Ltot,Lfib=Lfib,Lbase=Lbase,hstruct=hstruct,pi_j=pi_j,sub=sub,e_ref=e_ref)
}

set.seed(20)
A_true_A<-rnorm(15,0,1); A_true_A[c(1,2,13,15)]<-c(3.4,-3.1,3.0,-3.6)
resA<-run_scenario(c(1,1,1,2,2,3,4,6,10,15,25,2,1,8,1), A_true_A, "A: 15 groups (sigma_u well identified)")
resB<-run_scenario(c(1,1,1,2,2,3), c(3.8,-3.5,3.2,0.3,-0.4,0.2), "B: 6 groups (sigma_u barely identified)")

report<-function(r){ cat(sprintf("\n--- %s ---\n",r$label))
  cat(sprintf("  k_full: max=%.2f #(>0.7)=%d  ->  RB k_base: max=%.2f #(>0.7)=%d   (mean %.2f -> %.2f)\n",
    max(r$k_full),sum(r$k_full>0.7),max(r$k_base),sum(r$k_base>0.7),mean(r$k_full),mean(r$k_base)))
  cat(sprintf("  H1a k-hat predictor:  structural h_i=%.2f   influence L_i=%.2f  (Spearman vs k_full)\n", sp(r$hstruct,r$k_full), sp(r$Ltot,r$k_full)))
  cat(sprintf("  H2a base leverage -> base k-hat: Spearman(L_base,k_base)=%.2f\n", sp(r$Lbase,r$k_base)))
  cat(sprintf("  Acc vs refit (n=%d):  RB-LOO RMSE=%.3f   PSIS-LOO RMSE=%.3f\n",
    length(r$sub), sqrt(mean((r$e_rb[r$sub]-r$e_ref)^2)), sqrt(mean((r$e_full[r$sub]-r$e_ref)^2)))) }
cat("================ RESULTS ================"); report(resA); report(resB)

## ---- combined figure: rows = scenarios ----
png("leverage_psisloo.png", width=1550, height=820, res=135)
par(mfrow=c(2,4), mar=c(4,4.1,2.6,1))
panelset<-function(r){
  plot(r$hstruct,r$k_full,pch=19,col="#555555",xlab="structural leverage h_i (a-priori)",ylab="k-hat full",
       main=sprintf("k-hat predicted: rho=%.2f",sp(r$hstruct,r$k_full)),cex.main=.95); abline(h=0.7,lty=3,col="red")
  plot(r$k_full,r$k_base,pch=19,col=ifelse(r$Lbase>quantile(r$Lbase,.75),"#d62728","#2ca02c"),
       xlab="k-hat full (conditional)",ylab="k-hat base (RB)",main="RB-LOO cures fiber folds",cex.main=.95)
  abline(0,1,col="grey60",lty=2); abline(h=0.7,v=0.7,lty=3,col="grey70")
  plot(r$Lbase,r$k_base,pch=19,col="#d62728",xlab="base leverage L_base",ylab="k-hat base (RB)",
       main=sprintf("base leverage: rho=%.2f",sp(r$Lbase,r$k_base)),cex.main=.95); abline(h=0.7,lty=3,col="grey70")
  plot(r$e_ref,r$e_rb[r$sub],pch=19,col="#1f77b4",xlab="refit elpd (gold)",ylab="LOO elpd",main="RB (blue) vs PSIS (red) accuracy",cex.main=.95)
  points(r$e_ref,r$e_full[r$sub],pch=1,col="#d62728"); abline(0,1,lty=2,col="grey60")
}
panelset(resA); panelset(resB)
dev.off(); cat("\n\nWrote figure: leverage_psisloo.png  (top row = scenario A, bottom = B)\n")
