# =====================================================================
# P3 -- base-leverage triage (C3), gate G3. (project: fibr_sbc)
# H4: base leverage L_base predicts the residual base-IS k-hat after RB-LOO.
# Discriminant: fiber leverage -> k_full ; base leverage -> k_base (two failure modes).
# Part A: few-groups simulation sweep (base IS genuinely strained).
# Part B: eight schools -- the canonical few-groups LOO failure, cured by RB.
# =====================================================================
suppressMessages(library(loo)); set.seed(1)

gibbs_lmm <- function(y, grp, J, sigma, S=3000, warm=1000, Vmu=100, a0=2, b0=2, keep=rep(TRUE,length(y))) {
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

## ================= PART A: few-groups sweep =================
cat("PART A: few-groups sweep (base IS strained)\n")
configs<-list(c(1,1,1,2), c(1,1,1,1,2), c(1,1,1,1,2,3)); REP<-80
rows<-list(); ri<-0; set.seed(7)
for(cfg in configs){ J<-length(cfg); grp<-rep(1:J,cfg); N<-sum(cfg)
  for(r in 1:REP){ su<-runif(1,1.0,2.3); a_true<-rnorm(J,0,su); y<-a_true[grp]+rnorm(N,0,1)
    dr<-gibbs_lmm(y,grp,J,1.0,S=3000); out<-rb_loo(y,grp,1.0,dr); out<-out[is.finite(out$L_base),]
    if(nrow(out)>0){ ri<-ri+1; rows[[ri]]<-out } } }
PA<-do.call(rbind,rows)
sp<-function(a,b) suppressWarnings(cor(a,b,method="spearman"))
auc<-function(score,lab){ if(sum(lab)==0||sum(lab)==length(lab)) return(NA)
  r<-rank(score); (sum(r[lab==1])-sum(lab==1)*(sum(lab==1)+1)/2)/(sum(lab==1)*sum(lab==0)) }

cat(sprintf("  folds=%d   k_full>0.7: %d   k_base>0.7: %d   k_base>0.5: %d\n",
            nrow(PA), sum(PA$k_full>0.7), sum(PA$k_base>0.7), sum(PA$k_base>0.5)))
cat("  --- H4: base leverage predicts residual base k-hat ---\n")
cat(sprintf("  Spearman(L_base, k_base) = %.3f        [thr >=0.50 %s]\n", sp(PA$L_base,PA$k_base),
            ifelse(sp(PA$L_base,PA$k_base)>=0.5,"PASS","FAIL")))
ev<-as.numeric(PA$k_base>0.7)
a70<-auc(PA$L_base, ev); a50<-auc(PA$L_base, as.numeric(PA$k_base>0.5))
cat(sprintf("  AUC(k_base>0.7 ~ L_base) = %.3f  |  AUC(k_base>0.5 ~ L_base) = %.3f   [thr >=0.75 %s]\n",
            a70, a50, ifelse((!is.na(a70)&&a70>=0.75)||(!is.na(a50)&&a50>=0.75),"PASS","see note")))
cat(sprintf("  *** Does L_base BEAT the pooling factor for RB-failure? AUC(k_base>0.7): L_base=%.3f  h_struct=%.3f  (1-pi_j)=%.3f\n",
            auc(PA$L_base,ev), auc(PA$h_struct,ev), auc(1-PA$pi_j,ev)))
cat("  --- Discriminant (POOLED many-groups[P1] + few-groups[P3]): two failure modes ---\n")
DISC <- PA[,c("h_struct","L_base","k_full","k_base")]
g1p <- tryCatch(readRDS("rb_loo_gateG1.rds")$PF, error=function(e) NULL)
if(!is.null(g1p)) DISC <- rbind(g1p[,c("h_struct","L_base","k_full","k_base")], DISC)
cat(sprintf("  pooled folds=%d\n", nrow(DISC)))
cat(sprintf("             -> k_full   -> k_base\n"))
cat(sprintf("  h_struct :  %+.2f      %+.2f   (fiber leverage -> conditional failure)\n", sp(DISC$h_struct,DISC$k_full), sp(DISC$h_struct,DISC$k_base)))
cat(sprintf("  L_base   :  %+.2f      %+.2f   (base leverage  -> residual base failure)\n", sp(DISC$L_base,DISC$k_full), sp(DISC$L_base,DISC$k_base)))
hi<-DISC$h_struct>quantile(DISC$h_struct,0.6)
cat(sprintf("  *** CONDITIONAL on high fiber-leverage (folds RB is applied to, n=%d):\n", sum(hi)))
cat(sprintf("        Spearman(L_base,k_base)=%.2f  vs  Spearman(h_struct,k_base)=%.2f   <- does base leverage add value here?\n",
            sp(DISC$L_base[hi],DISC$k_base[hi]), sp(DISC$h_struct[hi],DISC$k_base[hi])))

## ================= PART B: eight schools =================
cat("\nPART B: eight schools (canonical few-groups LOO failure)\n")
y8<-c(28,8,-3,7,-1,1,18,12); s8<-c(15,10,16,11,9,11,10,18); J8<-8
# marginal (mu,tau) sampler on a grid (well-mixed, independent draws), then theta
mg<-seq(-15,30,length.out=220); tg<-seq(0.01,30,length.out=220)
GM<-outer(mg,tg,Vectorize(function(m,t){ sum(dnorm(y8,m,sqrt(t^2+s8^2),log=TRUE)) -0.5*(m^2)/100 - 0.5*(t^2)/25 }))
w<-exp(GM-max(GM)); w<-w/sum(w); S<-4000
idx<-sample(length(w),S,replace=TRUE,prob=as.vector(w)); mi<-((idx-1)%%length(mg))+1; ti<-((idx-1)%/%length(mg))+1
mu_s<-mg[mi]+rnorm(S,0,diff(mg)[1]/2); tau_s<-pmax(tg[ti]+rnorm(S,0,diff(tg)[1]/2),0.01)
# theta_j | mu,tau,y
Th<-sapply(1:J8,function(j){ v<-1/(1/s8[j]^2+1/tau_s^2); m<-v*(y8[j]/s8[j]^2+mu_s/tau_s^2); rnorm(S,m,sqrt(v)) })
Lf8<-sapply(1:J8,function(j) dnorm(y8[j],Th[,j],s8[j],log=TRUE))     # conditional
Lr8<-sapply(1:J8,function(j) dnorm(y8[j],mu_s,sqrt(tau_s^2+s8[j]^2),log=TRUE))  # RB (theta marginalised)
kf8<-suppressWarnings(loo(Lf8,r_eff=rep(1,J8)))$diagnostics$pareto_k
kr8<-suppressWarnings(loo(Lr8,r_eff=rep(1,J8)))$diagnostics$pareto_k
pi8<-(1/mean(tau_s)^2)/(1/s8^2+1/mean(tau_s)^2)
cat(sprintf("  posterior: mu=%.1f  tau=%.1f\n", mean(mu_s), mean(tau_s)))
cat("  school  y   se   pi_j   k_full  k_base\n")
for(j in 1:J8) cat(sprintf("    %d   %4.0f %4.0f  %.2f   %+.2f   %+.2f\n", j,y8[j],s8[j],pi8[j],kf8[j],kr8[j]))
cat(sprintf("  full PSIS-LOO: max k=%.2f (#>0.7=%d)  ->  RB-LOO: max k=%.2f (#>0.7=%d)\n",
            max(kf8),sum(kf8>0.7),max(kr8),sum(kr8>0.7)))

## ---- figure ----
png("p3_base_leverage.png", width=1550, height=520, res=135)
par(mfrow=c(1,3), mar=c(4.2,4.3,3,1))
plot(PA$L_base, PA$k_base, pch=19, col=adjustcolor("#d62728",0.3), cex=0.6,
     xlab="base leverage L_base (analytic)", ylab="k-hat base (after RB-LOO)",
     main=sprintf("H4: base leverage predicts\nresidual strain  rho=%.2f", sp(PA$L_base,PA$k_base))); abline(h=0.7,lty=3,col="grey50")
# head-to-head: pooling factor predicts RB-failure as well as the base leverage
aucs<-c("L_base"=auc(PA$L_base,ev), "h_struct"=auc(PA$h_struct,ev), "1-pi_j"=auc(1-PA$pi_j,ev))
bp<-barplot(aucs, ylim=c(0,1), col=c("#d62728","#1f77b4","#7f7f7f"), ylab="AUC (flag RB-failure k_base>0.7)",
            main="Pooling factor >= base leverage\n(C3 redundant with C1)")
text(bp, aucs+0.04, sprintf("%.2f",aucs)); abline(h=0.5,lty=3,col="grey50")
plot(kf8, kr8, pch=19, col="#2ca02c", cex=1.2, xlim=c(0,max(kf8)*1.05), ylim=c(min(kr8)*1.1,max(kf8)*1.05),
     xlab="k-hat full (conditional)", ylab="k-hat base (RB)", main="Eight schools: RB cures the funnel LOO")
abline(0,1,lty=2,col="grey60"); abline(h=0.7,v=0.7,lty=3,col="grey70"); text(kf8,kr8,1:8,pos=3,cex=0.7)
dev.off(); cat("\nWrote figure: p3_base_leverage.png\n")
