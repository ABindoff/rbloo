# =====================================================================
# RB-PIT in the bjlm setting: is the analytic conditional calibrated for the
# NON-LINEAR smooth change-point coordinate, and where does it break? (fibr_sbc)
#
# mu(tau) = b0 + b1*(tau-omega) + delta*(tau-omega)*sigma(rho*(tau-omega))   [bjlm]
# Random change-point omega ~ N(omega0, sig_om^2); Gaussian outcome, sigma known.
# pooling_factor.R linearises via the analytic gradient
#   dmu/domega = -( b1 + delta*s*(1 + d*rho*(1-s)) ),  d=tau-omega, s=sigma(rho*d)
# -> Fisher precision Glik = sum(dmu^2)/sigma^2.
#
# Controlled isolation (all other params known = best case for RB). Per draw:
#   omega~ ~ prior;  y_i = mu(tau_i;omega~)+N(0,sigma^2);  exact conditional p(omega|y)
#   on a dense grid.  Three PITs vs uniform:
#     - exact          : positive control, uniform by construction
#     - moment-matched : best possible Gaussian (exact cond. mean+sd); its
#                        miscalibration = pure NON-GAUSSIANITY of the conditional
#     - Fisher / pf    : bjlm's approximation, mode + Fisher precision (Glik)
# Sweeps: identification strength (n_obs) and change-point sharpness (rho).
# =====================================================================
set.seed(7)
sigm  <- function(x) 1/(1+exp(-x))
mu_cp <- function(tau, omega, b0, b1, delta, rho){d<-tau-omega; b0+b1*d+delta*d*sigm(rho*d)}
dmu_domega <- function(tau, omega, b1, delta, rho){d<-tau-omega; s<-sigm(rho*d); -(b1+delta*s*(1+d*rho*(1-s)))}
b0<-0; b1<--0.3; delta<-1.2; sigma<-0.4; omega0<-3; sig_om<-0.8

analyse_one <- function(omega_true, rho, taus){
  y  <- mu_cp(taus, omega_true, b0,b1,delta,rho) + rnorm(length(taus),0,sigma)
  gr <- seq(omega0-6*sig_om, omega0+6*sig_om, length.out=2000); dx<-gr[2]-gr[1]
  lp <- dnorm(gr,omega0,sig_om,log=TRUE) +
        vapply(gr,function(w) sum(dnorm(y,mu_cp(taus,w,b0,b1,delta,rho),sigma,log=TRUE)),numeric(1))
  dens <- exp(lp-max(lp)); dens <- dens/sum(dens*dx); cdf<-cumsum(dens)*dx; cdf<-cdf/max(cdf)
  pit_exact <- approx(gr,cdf,omega_true,rule=2)$y
  m_mm <- sum(gr*dens*dx); s_mm <- sqrt(sum((gr-m_mm)^2*dens*dx))       # exact cond moments
  pit_mm <- pnorm(omega_true, m_mm, s_mm)
  mode <- gr[which.max(lp)]                                             # Fisher/pf: mode + Glik
  Glik <- sum(dmu_domega(taus,mode,b1,delta,rho)^2)/sigma^2
  pit_fis <- pnorm(omega_true, mode, 1/sqrt(1/sig_om^2+Glik))
  c(exact=pit_exact, mm=pit_mm, fisher=pit_fis)
}
ks <- function(u){u<-sort(u);N<-length(u);i<-1:N;max(pmax(i/N-u,u-(i-1)/N))}
ks_crit <- function(N,M=3000){D<-numeric(M);for(m in 1:M){u<-sort(runif(N));i<-1:N;D[m]<-max(pmax(i/N-u,u-(i-1)/N))};quantile(D,.95)}

S<-800; crit<-as.numeric(ks_crit(S))
sweep_run <- function(rho, n_obs){taus<-seq(0.3,5.7,length.out=n_obs)
  t(vapply(1:S, function(s) analyse_one(rnorm(1,omega0,sig_om),rho,taus), numeric(3)))}

cat(sprintf("bjlm change-point RB, S=%d, 0.95 KS crit=%.4f\n\n",S,crit))
n_vec <- c(3,4,6,8,12); rho_fix <- 4
cat("Sweep A: identification strength (n_obs), rho=4\n")
resN <- lapply(n_vec, function(n) sweep_run(rho_fix,n)); names(resN)<-n_vec
for(i in seq_along(n_vec)) cat(sprintf("  n_obs=%2d: KS exact=%.4f  moment-match=%.4f  Fisher/pf=%.4f  %s\n",
  n_vec[i],ks(resN[[i]][,1]),ks(resN[[i]][,2]),ks(resN[[i]][,3]),
  ifelse(ks(resN[[i]][,3])>crit,"<- Fisher/pf OUT","in band")))
rho_vec <- c(1,3,8,20); n_fix<-8
cat("Sweep B: sharpness (rho), n_obs=8\n")
resR <- lapply(rho_vec, function(r) sweep_run(r,n_fix)); names(resR)<-rho_vec
for(i in seq_along(rho_vec)) cat(sprintf("  rho=%2d: KS exact=%.4f  moment-match=%.4f  Fisher/pf=%.4f  %s\n",
  rho_vec[i],ks(resR[[i]][,1]),ks(resR[[i]][,2]),ks(resR[[i]][,3]),
  ifelse(ks(resR[[i]][,3])>crit,"<- Fisher/pf OUT","in band")))

# ------------------------------ figure ------------------------------
ed<-function(u){u<-sort(u);list(x=u,y=(1:length(u))/length(u)-u)}
panel<-function(P,main){yl<-max(crit*1.3,max(abs(ed(P[,2])$y)),max(abs(ed(P[,3])$y)))
  plot(NA,xlim=c(0,1),ylim=c(-yl,yl),xlab="PIT",ylab=expression(hat(F)-u),main=main,cex.main=0.94)
  polygon(c(0,1,1,0),c(crit,crit,-crit,-crit),col=adjustcolor("grey50",.13),border=NA);abline(h=0,col="grey70")
  e<-ed(P[,1]);lines(e$x,e$y,type="s",col="grey55",lwd=1.5)
  e<-ed(P[,2]);lines(e$x,e$y,type="s",col="#2ca02c",lwd=1.5)
  e<-ed(P[,3]);lines(e$x,e$y,type="s",col="#d62728",lwd=1.9)
  legend("topleft",bty="n",cex=0.78,col=c("grey55","#2ca02c","#d62728"),lwd=2,legend=c(
    sprintf("exact ctrl D=%.3f",ks(P[,1])),sprintf("moment-match D=%.3f",ks(P[,2])),sprintf("Fisher/pf D=%.3f",ks(P[,3]))))}
png("bjlm_changepoint_rb.png",width=1500,height=1000,res=140)
par(mfrow=c(2,3),mar=c(4,4.2,3,1))
tg<-seq(0,6,length.out=200); cols<-c("#1f77b4","#ff7f0e","#d62728")
plot(NA,xlim=c(0,6),ylim=range(sapply(c(1,4,20),function(r)mu_cp(tg,omega0,b0,b1,delta,r))),xlab=expression(tau),ylab=expression(mu(tau)),main="bjlm mean sharpens with rho")
for(i in seq_along(c(1,4,20))) lines(tg,mu_cp(tg,omega0,b0,b1,delta,c(1,4,20)[i]),col=cols[i],lwd=2)
abline(v=omega0,lty=3,col="grey50");legend("topright",bty="n",legend=paste0("rho=",c(1,4,20)),col=cols,lwd=2)
panel(resN[["3"]],"Sparse: n_obs=3 (weak id)")
panel(resN[["8"]],"Well-identified: n_obs=8")
# KS vs n_obs
kN<-sapply(resN,function(P)c(ks(P[,2]),ks(P[,3])))
plot(n_vec,kN[2,],type="b",pch=19,col="#d62728",ylim=c(0,max(kN,crit)*1.1),lwd=2,xlab="observations n_obs",ylab="PIT miscalibration (KS D)",main="RB error vs identification")
lines(n_vec,kN[1,],type="b",pch=17,col="#2ca02c",lwd=2);abline(h=crit,lty=2,col="grey50")
legend("topright",bty="n",cex=0.85,legend=c("Fisher/pf","moment-match","0.95 band"),col=c("#d62728","#2ca02c","grey50"),lwd=2,pch=c(19,17,NA),lty=c(1,1,2))
# KS vs rho
kR<-sapply(resR,function(P)c(ks(P[,2]),ks(P[,3])))
plot(rho_vec,kR[2,],type="b",pch=19,col="#d62728",ylim=c(0,max(kR,crit)*1.1),lwd=2,log="x",xlab="sharpness rho (log)",ylab="PIT miscalibration (KS D)",main="Sharpness alone: RB stays calibrated")
lines(rho_vec,kR[1,],type="b",pch=17,col="#2ca02c",lwd=2);abline(h=crit,lty=2,col="grey50")
legend("topleft",bty="n",cex=0.85,legend=c("Fisher/pf","moment-match","0.95 band"),col=c("#d62728","#2ca02c","grey50"),lwd=2,pch=c(19,17,NA),lty=c(1,1,2))
dev.off();cat("\nWrote figure: bjlm_changepoint_rb.png\n")
