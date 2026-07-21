# =====================================================================
# Control-variate PIT for SBC  --  Phase 1b prototype  (project: fibr_sbc)
#
# Question: can we keep the sampler's fiber draws (so we still TEST them, unlike
# replace-RB which is blind) AND cut variance toward replace-RB, using the fibr
# conditional CDF as a control variate?
#
#   u_cv = (1/M) sum_{m=1..M} c(base_m)  +  (1/L) sum_{l in sub} ( g_l - c_l )
#   g_l = 1{ alpha_j^(l) < alpha_j~ }     c = Phi( (alpha_j~ - m_j)/s_j )   [fibr]
#
# with M>>L emulating the target model's cost asymmetry (c cheap on the low-dim
# base; fiber = latent GP field, expensive).
#
# ANSWER (this prototype): NO free lunch. At small L the dominant Monte-Carlo
# noise is the fiber-conditional SAMPLING noise Var(g|base)=c(1-c), which the CV
# KEEPS (that is what makes it fiber-sensitive). It removes only the base-marginal
# component Var(c), small when the base is well identified. So u_cv variance ~=
# counting, NOT ~= replace-RB. Replace-RB's low variance comes precisely from
# discarding g -> blindness. You cannot have both at small L.
#
# What CV *does* buy honestly: it is CONTINUOUS (no rank discreteness) and it does
# retain fiber-bug sensitivity. Useful, but it is not the base-side L-slashing win.
# =====================================================================
set.seed(11)
fibr_moments <- function(pp, li){prec<-pp+li; list(pi=pp/prec, s=1/sqrt(prec))}
J<-6; n_j<-c(2,5,10,20,40,80); sigma_y<-1; tau<-5; a0<-3; b0<-2

sim_prior <- function(){mu<-rnorm(1,0,tau); s2<-1/rgamma(1,a0,b0)
  al<-rnorm(J,mu,sqrt(s2)); yb<-rnorm(J,al,sigma_y/sqrt(n_j)); list(mu=mu,sigma2=s2,alpha=al,ybar=yb)}
gibbs <- function(d,M=800,warmup=200,alpha_infl=1,s2_bias=1){yb<-d$ybar;mu<-0;s2<-1
  km<-numeric(M);ks<-numeric(M);ka<-matrix(0,M,J);tot<-warmup+M
  for(t in 1:tot){pa<-1/s2+n_j/sigma_y^2;va<-1/pa;ma<-va*(mu/s2+n_j*yb/sigma_y^2)
    al<-rnorm(J,ma,sqrt(va)*alpha_infl);vm<-1/(J/s2+1/tau^2);mu<-rnorm(1,vm*sum(al)/s2,sqrt(vm))
    s2<-s2_bias/rgamma(1,a0+J/2,b0+0.5*sum((al-mu)^2))
    if(t>warmup){i<-t-warmup;km[i]<-mu;ks[i]<-s2;ka[i,]<-al}}
  list(mu=km,sigma2=ks,alpha=ka,M=M)}

pits <- function(ch,d,L){M<-ch$M;idx<-round(seq(1,M,length.out=L))
  o<-matrix(0,J,3,dimnames=list(NULL,c("cnt","rb","cv")));pp<-1/ch$sigma2
  for(j in 1:J){mo<-fibr_moments(pp,n_j[j]/sigma_y^2);m<-mo$pi*ch$mu+(1-mo$pi)*d$ybar[j]
    ca<-pnorm((d$alpha[j]-m)/mo$s);g<-as.numeric(ch$alpha[idx,j]<d$alpha[j]);cs<-ca[idx]
    o[j,]<-c(mean(g),mean(ca),min(max(mean(ca)+mean(g-cs),0),1))}
  o}
ks_stat<-function(u){u<-sort(u);N<-length(u);i<-1:N;max(pmax(i/N-u,u-(i-1)/N))}
ks_crit<-function(N,l=.05,Mc=3000){D<-numeric(Mc);for(m in 1:Mc){u<-sort(runif(N));i<-1:N;D[m]<-max(pmax(i/N-u,u-(i-1)/N))};quantile(D,1-l)}
ks_cdisc<-function(N,Ld,l=.05,Mc=3000){D<-numeric(Mc);for(m in 1:Mc){u<-sample(0:Ld,N,replace=TRUE)/Ld;D[m]<-ks_stat(u)};quantile(D,1-l)}
run<-function(S,L,M=800,ai=1,sb=1){A<-array(0,c(S,J,3));for(s in 1:S){d<-sim_prior();ch<-gibbs(d,M=M,alpha_infl=ai,s2_bias=sb);A[s,,]<-pits(ch,d,L)};list(cnt=as.vector(A[,,1]),rb=as.vector(A[,,2]),cv=as.vector(A[,,3]))}

L<-25; Mb<-800
cat(sprintf("Cost asymmetry: M=%d cheap base draws, L=%d expensive fiber draws\n\n",Mb,L))
cat("E1 validity (correct sampler), S=400\n")
e1<-run(400,L,Mb); cc<-as.numeric(ks_crit(length(e1$cnt))); cd<-as.numeric(ks_cdisc(length(e1$cnt),L))
cat(sprintf("  D: cnt=%.4f (disc crit %.4f) | rb=%.4f | cv=%.4f (cont crit %.4f)%s\n",
  ks_stat(e1$cnt),cd,ks_stat(e1$rb),ks_stat(e1$cv),cc,
  ifelse(ks_stat(e1$cv)>cc,"  <- CV slightly over-dispersed (within-sim noise)","")))

cat("E2 within-sim variance (fixed dataset, 200 reruns)\n")
d0<-sim_prior(); R<-200; Mc_<-matrix(0,R,J);Mr<-matrix(0,R,J);Mv<-matrix(0,R,J)
for(r in 1:R){ch<-gibbs(d0,M=Mb);p<-pits(ch,d0,L);Mc_[r,]<-p[,1];Mr[r,]<-p[,2];Mv[r,]<-p[,3]}
sd_c<-apply(Mc_,2,sd);sd_r<-apply(Mr,2,sd);sd_v<-apply(Mv,2,sd)
for(j in 1:J)cat(sprintf("  n=%2d: SD cnt=%.4f rb=%.4f cv=%.4f  (cv reduces cnt var x%.2f; rb x%.0f)\n",
  n_j[j],sd_c[j],sd_r[j],sd_v[j],(sd_c[j]/sd_v[j])^2,(sd_c[j]/sd_r[j])^2))

cat("E3 fiber bug (alpha infl=0.8), S=400\n")
e3<-run(400,L,Mb,ai=0.8)
cat(sprintf("  D: cnt=%.4f | rb=%.4f (BLIND) | cv=%.4f\n",ks_stat(e3$cnt),ks_stat(e3$rb),ks_stat(e3$cv)))
cat("E4 base bug (s2 bias=0.7), S=400\n")
e4<-run(400,L,Mb,sb=0.7)
cat(sprintf("  D: cnt=%.4f | rb=%.4f | cv=%.4f\n",ks_stat(e4$cnt),ks_stat(e4$rb),ks_stat(e4$cv)))

# variance decomposition on the fixed dataset (WHY CV cannot win)
cat("E5 variance decomposition (fixed dataset, long chain)\n")
chd<-gibbs(d0,M=4000,warmup=500); pp<-1/chd$sigma2
Vbet<-Ewi<-numeric(J)
for(j in 1:J){mo<-fibr_moments(pp,n_j[j]/sigma_y^2);m<-mo$pi*chd$mu+(1-mo$pi)*d0$ybar[j]
  ca<-pnorm((d0$alpha[j]-m)/mo$s);Vbet[j]<-var(ca);Ewi[j]<-mean(ca*(1-ca))}
for(j in 1:J)cat(sprintf("  n=%2d: base-marginal Var(c)=%.4f (CV removes) | fiber-noise E[c(1-c)]=%.4f (CV keeps)\n",
  n_j[j],Vbet[j],Ewi[j]))

# =============================== FIGURE ===============================
ed<-function(u){u<-sort(u);list(x=u,y=(1:length(u))/length(u)-u)}
col_c<-"#1f77b4";col_r<-"#7f7f7f";col_v<-"#d62728"
overlay<-function(lst,crit,main){yl<-max(crit*1.4,max(sapply(lst,function(u)max(abs(ed(u)$y)))))
  plot(NA,xlim=c(0,1),ylim=c(-yl,yl),xlab="PIT value",ylab=expression(hat(F)(u)-u),main=main,cex.main=0.95)
  polygon(c(0,1,1,0),c(crit,crit,-crit,-crit),col=adjustcolor("grey50",.13),border=NA);abline(h=0,col="grey70")
  cols<-c(col_c,col_r,col_v);for(i in seq_along(lst)){e<-ed(lst[[i]]);lines(e$x,e$y,type="s",col=cols[i],lwd=1.7)}
  legend("topleft",bty="n",cex=0.8,col=cols,lwd=2,legend=c(
    sprintf("counting  D=%.3f",ks_stat(lst[[1]])),sprintf("replace-RB D=%.3f",ks_stat(lst[[2]])),
    sprintf("CV-PIT  D=%.3f",ks_stat(lst[[3]]))))}
png("cv_pit.png",width=1500,height=1450,res=140)
par(mfrow=c(3,2),mar=c(4,4.3,3,1))
overlay(list(e1$cnt,e1$rb,e1$cv),cc,"E1  Correct sampler")
bp<-barplot(rbind(sd_c,sd_r,sd_v),beside=TRUE,col=c(col_c,col_r,col_v),names.arg=paste0("n=",n_j),
  ylab="within-sim SD of PIT",main="E2  CV variance ~= counting, NOT replace-RB",cex.names=.8,ylim=c(0,0.13))
legend("top",bty="n",horiz=TRUE,cex=.8,legend=c("counting","replace-RB","CV"),fill=c(col_c,col_r,col_v))
overlay(list(e3$cnt,e3$rb,e3$cv),cc,"E3  Fiber bug: RB BLIND, CV & counting detect")
overlay(list(e4$cnt,e4$rb,e4$cv),cc,expression("E4  Base ("*sigma^2*") bug: all detect"))
barplot(rbind(Ewi,Vbet),beside=FALSE,col=c(col_v,col_c),names.arg=paste0("n=",n_j),
  ylab="per-draw variance",main="E5  Why: fiber noise (red) dominates, CV can't cut it",cex.names=.8)
legend("topright",bty="n",cex=.8,legend=c("fiber-sampling noise E[c(1-c)] (CV keeps)","base-marginal Var(c) (CV removes)"),fill=c(col_v,col_c))
plot.new();legend("topleft",bty="n",cex=.92,title="Verdict",legend=c(
 "No free lunch at small L:",
 "* replace-RB: min variance, BLIND to fiber (E2,E3)",
 "* counting: tests fiber, noisy + discrete",
 "* CV-PIT: tests fiber (E3), continuous, but variance",
 "   ~= counting (E2) -- it keeps the fiber Bernoulli",
 "   noise E[c(1-c)] that dominates at small L (E5).",
 "",
 "Design: base coords -> replace-RB (slash L).",
 "Fiber coords -> need real fiber draws; CV adds",
 "continuity + removes base noise when base diffuse."))
dev.off();cat("\nWrote figure: cv_pit.png\n")
