
#########################################################################
#									
# 	This R Program simulates the method of Wages and Tait (2015)* 				
#									
#	*Wages, N. and Tait, C. (2015) Seamless Phase I/II Adaptive Design 
#	for Oncology Trials of Molecularly Targeted Therapies. J Biopharm Stats.							
#									
#########################################################################


###install required R packages
library(nnet)
library(dfcrm)
library(binom)


###Load the function 'bpocrm' 
bpocrm<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb){
  
  
  # if a single ordering is inputed as a vector, convert it to a matrix
  if(is.vector(q.skel)) q.skel=t(as.matrix(q.skel));
  
  nord.eff = nrow(q.skel);
  mprior.eff = rep(1/nord.eff, nord.eff); # prior for each efficacy ordering
  
  avar=1.34
  bcrmh<-function(a,p,y,n){
    s2=avar
    lik=exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  bcrmht<-function(a,p,y,n){
    s2=avar
    lik=a*exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  bcrmht2<-function(a,p,y,n){
    s2=avar
    lik=a^2*exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  
  ### run a trial 	
  ncomb = ncol(q.skel);   #number of combos
  y=rep(0,ncomb);  #number of toxicity/responses at each dose level
  z=rep(0,ncomb);   #number of efficacy at each dose level
  n=rep(0,ncomb);  #number of treated patients at each dose level
  comb.curr = start.comb;  # current dose level	 
  ptox.hat = numeric(ncomb); # estimate of toxicity prob
  comb.select=rep(0,ncomb); # a vector of indicators for dose selection
  stop=stops=stopf=0; #indicate if trial stops early
  i=1	
  while(i <= ncohort)
  {
    # generate data for a new cohort of patients
    
    y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,p0[comb.curr]);
    z[comb.curr] = z[comb.curr] + rbinom(1,cohortsize,q0[comb.curr]);
    n[comb.curr] = n[comb.curr] + cohortsize;		
    
    ####Toxicity
    marginal= integrate(bcrmh,lower=-Inf,upper=Inf, p=p.skel, y=y,n=n,abs.tol = 0)$value;
    est=integrate(bcrmht,lower=-10,upper=10, p.skel, y, n,abs.tol = 0)$value/marginal
    ptox.hat=p.skel**exp(est)
    
    
    ####Efficacy
    marginal.eff = est.eff=rep(0, nord.eff);
    for(k in 1:nord.eff)
    {
      marginal.eff[k] = integrate(bcrmh,lower=-Inf,upper=Inf, p=q.skel[k,], y=z,n=n,abs.tol = 0)$value;
      est.eff[k]=integrate(bcrmht,lower=-10,upper=10, q.skel[k,], z, n,abs.tol = 0)$value/marginal.eff[k]
    }		
    postprob.eff = (marginal.eff*mprior.eff)/sum(marginal.eff*mprior.eff);
    # efficacy model selection, identify the model with the highest posterior prob
    if(nord.eff>1){ 
      meff.sel = which.is.max(postprob.eff); 
    } else{
      meff.sel = 1;
    }
    
    aset=which(ptox.hat<=tul)
    if(length(aset)==0){aset=which.min(ptox.hat)}
    peff.hat=matrix(,nrow=nrow(q.skel),ncol=d)
    peff.hat=q.skel**exp(est.eff)
    
    peff.hat.aset=rep(0,ncomb)
    peff.hat.aset[aset]=peff.hat[meff.sel,aset]
    ri0=peff.hat.aset/sum(peff.hat.aset)
    
    if(length(aset)==1){
      comb.best=aset
    } else {
      ifelse(sum(n)<=n.ar,comb.best<-sample(1:ncomb,1,prob=ri0),comb.best<-which.max(peff.hat.aset))
    } 		
    
    ##########stopping rules
    safety=binom.confint(y[1],n[1],conf.level=0.95,methods="exact")$lower
    if(safety>tul){
      stops=1
      break
    }
    
    if(sum(n) > n.ar){
      futility=binom.confint(z[comb.best],n[comb.best],conf.level=0.95,methods="exact")$upper
      if(futility<ell){
        stopf=1
        break
      }
    } 
    
    comb.curr<-comb.best
    i=i+1
  }
  
  if(stop==0 & stops==0 & stopf==0){
    comb.select[comb.curr]=comb.select[comb.curr]+1;
  }
  return(list(comb.select=comb.select,tox.data=y,eff.data=z,pt.allocation=n,stopf=stopf,stops=stops))
}
##########'bpocrm' end here

###Load the function 'bpocrm.sim' 
bpocrm.sim<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb){
  ncomb=length(p0)
  
  comb.select<-y<-z<-n<-matrix(nrow=ntrial,ncol=ncomb)
  nstopf=nstops=0
  
  for(i in 1:ntrial){
    result<-bpocrm(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb)
    comb.select[i,]=result$comb.select
    y[i,]=result$tox.data
    z[i,]=result$eff.data
    n[i,]=result$pt.allocation
    nstopf=nstopf+result$stopf
    nstops=nstops+result$stops
  }
  cat("True tox probability:           ", round(p0,3), sep="\t",  "\n");
  cat("True eff probability:           ", round(q0,3), sep="\t",  "\n");
  cat("selection percentage:           ", formatC(colMeans(comb.select)*100, digits=1, format="f"), sep="\t",  "\n");
  cat("number of toxicities:           ", formatC(colMeans(y), digits=1, format="f"), sep="\t",   "\n");
  cat("number of responses:            ", formatC(colMeans(z), digits=1, format="f"), sep="\t",   "\n");
  cat("number of patients treated:     ", formatC(colMeans(n), digits=1, format="f"), sep="\t",   "\n");
  cat("percentage of stop (safety):    ", nstops/ntrial*100, "\n");
  cat("percentage of stop (futility):  ", nstopf/ntrial*100, "\n");
}
##########'bpocrm.sim' end here


######################################
#
#        Generate results in Table 3
#              for psi=0
#
######################################

#####Specify the total number of doses
d<-5

###Specify a set of toxicity skeleton values
p.skel<-c(0.01,0.08,0.15,0.22,0.29)

#####Specify the number of possible efficacy orderings
g<-9   #efficacy

###Specifiy the possible efficacy orderings of the doses
q.skel<-matrix(nrow=g,ncol=d)
q.skel[1,]<-c(0.60,0.70,0.60,0.50,0.40)
q.skel[2,]<-c(0.70,0.60,0.50,0.40,0.30)
q.skel[3,]<-c(0.50,0.60,0.70,0.60,0.50)
q.skel[4,]<-c(0.40,0.50,0.60,0.70,0.60)
q.skel[5,]<-c(0.30,0.40,0.50,0.60,0.70)
q.skel[6,]<-c(0.70,0.70,0.70,0.70,0.70) 
q.skel[7,]<-c(0.60,0.70,0.70,0.70,0.70)
q.skel[8,]<-c(0.50,0.60,0.70,0.70,0.70) 
q.skel[9,]<-c(0.40,0.50,0.60,0.70,0.70)

p1<-c(0.01,0.05,0.10,0.15,0.20)
q1<-c(0.30,0.50,0.60,0.40,0.25)

p2<-c(0.02,0.06,0.12,0.30,0.40)
q2<-c(0.38,0.50,0.40,0.30,0.25)

p3<-c(0.03,0.09,0.16,0.28,0.42)
q3<-c(0.25,0.35,0.48,0.65,0.52)

p4<-c(0.02,0.05,0.07,0.09,0.11)
q4<-c(0.68,0.56,0.49,0.40,0.33)



tul<-0.33    ##toxicity upper limit 
ell<-0.20    ##efficacy lower limit
cohortsize=1 ##cohort size for each inclusion
ncohort=48   ##number of cohorts
start.comb=1 ##starting dose
n.ar=24      ##size of AR phase
ntrial=500   ##number of simulated trials 

p0<-p4
q0<-q4
set.seed(580)    ##random seed
##simulate one trial
#bpocrm(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb)

##simulate many trials
bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb)





