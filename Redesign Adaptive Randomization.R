
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

########################Functions in the BinarySimCLF package (the package is no longer available on CRAN)##########################################################
`ranBin2` <-
  function(nRep,u,psi,seed)
  {
    if (!missing(seed))
    {
      set.seed(seed);
    }
    u12 <- solve2(u[1], u[2], psi);
    y <- matrix( rep(-1,2*nRep),nrow=nRep );
    y[,1] <- ifelse( runif(nRep)<=u[1], 1, 0 );
    y[,2] <- y[,1]*(runif(nRep) <= u12/u[1]) + (1-y[,1])*(runif(nRep) <= (u[2]-u12)/(1-u[1]));
    return(y)
  }
`solve2` <-
  function(mui, muj, psi)
  {
    if (psi == 1)
    {
      return(mui*muj);
    }
    else if (psi != 1)
    {
      a <- 1 - psi;
      b <- 1 - a*(mui + muj);
      c <- -psi*(mui * muj);
      muij <- mardia(a,b,c);
    }
    return(muij);
  }
`mardia` <-
  function(a,b,c)
  {
    if (a == 0 )
    {
      return(-c/b);
      
    }
    k <- ifelse( b > 0, 1, ifelse(b < 0, -1, 0) );
    p <- -0.5 * (b + k*sqrt(b^2-4*a*c));
    r1 <- p/a;
    r2 <- c/p;
    
    r <- ifelse( r2 > 0, r2, r1 );
    return(r);
  }

#End of binarySinCLF#################################################################################
###Load the function 'bpocrm' 
bpocrm<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb,ar.strategy){
  
  
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
    
    #y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,p0[comb.curr]);
    #z[comb.curr] = z[comb.curr] + rbinom(1,cohortsize,q0[comb.curr]);
    #n[comb.curr] = n[comb.curr] + cohortsize;		
    #Correlated toxicity and efficacy #################################################
    joint.prob <- c(p0[comb.curr],q0[comb.curr]) 
    joint.outcome <- ranBin2(cohortsize,joint.prob,psi)
    y[comb.curr] = y[comb.curr] + joint.outcome[,1];
    z[comb.curr] = z[comb.curr] + joint.outcome[,2];
    n[comb.curr] = n[comb.curr] + cohortsize
    ###################################################################################
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
    ######################################################Redesign adaptive randomization##########################
    #Original: 
    if (ar.strategy==0)
    {
      ar.prob=peff.hat.aset/sum(peff.hat.aset)
      
     
    if(length(aset)==1){
      comb.best=aset
    } else {
      ifelse(sum(n)<=n.ar,comb.best<-sample(1:ncomb,1,prob=ar.prob),comb.best<-which.max(peff.hat.aset))
    }
      
    }#end of original strategy
    ##################################################################################################################
    
    #Alternative Randomization Strategy 1:
    if (ar.strategy==1) {
      ar.prob <- rep(0,ncomb)
      for (j in 1:length(postprob.eff)){
        skel.best <- which.max(q.skel[j,]) #The best dose by each model
        if (!is.element(skel.best,aset)){skel.best <- max(aset)}
        skel.best <- ifelse(is.element(skel.best,aset),skel.best,max(aset))
        ar.prob[skel.best] <- ar.prob[skel.best]+postprob.eff[j]
      } 
      if(length(aset)==1){
        comb.best=aset
      } else {
        ifelse(sum(n)<=n.ar,comb.best<-sample(1:ncomb,1,prob=ar.prob),comb.best<-which.max(peff.hat.aset))
      }
    }#End of strategy 1
    
    # ALternative Randomization Strategy 2:
    if (ar.strategy==2){
      ar.prob <- rep(0,ncomb)
      for (j in 1:length(postprob.eff)){
        skel.best <- which.max(q.skel[j,]) #The best dose by each model
        
        skel.best <- ifelse(is.element(skel.best,aset),skel.best,max(aset))
        
        if (postprob.eff[j] >= 1/length(postprob.eff)){
          ar.prob[skel.best] <- ar.prob[skel.best]+postprob.eff[j]
        } else {
          ar.prob[skel.best] <- ar.prob[skel.best]+0
        }
      } 
      ar.prob <- ar.prob/sum(ar.prob) 
      if(length(aset)==1){
        comb.best=aset
      } else {
        ifelse(sum(n)<=(cohortsize*ncohort-cohortsize),comb.best<-sample(1:ncomb,1,prob=ar.prob),comb.best<-which.max(peff.hat.aset))
        #unless it is the last sample(cohort), we keep randomizing
      }
    }
    ######################################
    # END of strategy 2
    # ALternative Randomization Strategy 3:
    if (ar.strategy==3){
      ar.prob <- rep(0,ncomb)
      #nord.eff is the number of candidate efficacy models
      nord.eff.limit<-ceiling(nord.eff*( (ncohort*cohortsize-sum(n))/ncohort*cohortsize))
      #nord.eff.limit is the number of models we will consider
      if (nord.eff.limit ==0) {
        cut.off=max(postprob.eff)
      } else {
        cut.off.order <- order(postprob.eff,decreasing=F)[nord.eff-nord.eff.limit+1]
        cut.off <- postprob.eff[cut.off.order]
      }
      
      for (j in 1:length(postprob.eff)){
        skel.best <- which.max(q.skel[j,]) #The best dose by each model
        if (!is.element(skel.best,aset)){skel.best <- max(aset)}
        skel.best <- ifelse(is.element(skel.best,aset),skel.best,max(aset))
        
        if (postprob.eff[j] >= cut.off){
          ar.prob[skel.best] <- ar.prob[skel.best]+postprob.eff[j]
        } else {
          ar.prob[skel.best] <- ar.prob[skel.best]+0
        }
      }
      ar.prob <- ar.prob/sum(ar.prob)    
      if(length(aset)==1){
        comb.best=aset
      } else {
        ifelse(sum(n)<=(cohortsize*ncohort-cohortsize),comb.best<-sample(1:ncomb,1,prob=ar.prob),comb.best<-which.max(peff.hat.aset))
        #unless it is the last sample(cohort), we keep randomizing
      }

    } #END of Strategy 3
    
    # ALternative Randomization Strategy 4:
    if (ar.strategy==4){
      ar.prob <- rep(0,ncomb)
      #nord.eff is the number of candidate efficacy models
      nord.eff.limit<-ceiling(nord.eff*(( (ncohort*cohortsize-sum(n))/ncohort*cohortsize))**drop.rate)
      #nord.eff.limit is the number of models we will consider
      if (nord.eff.limit ==0) {
        cut.off=max(postprob.eff)
      } else {
        cut.off.order <- order(postprob.eff,decreasing=F)[nord.eff-nord.eff.limit+1]
        cut.off <- postprob.eff[cut.off.order]
      }
      
      for (j in 1:length(postprob.eff)){
        skel.best <- which.max(q.skel[j,]) #The best dose by each model
        if (!is.element(skel.best,aset)){skel.best <- max(aset)}
        skel.best <- ifelse(is.element(skel.best,aset),skel.best,max(aset))
        
        if (postprob.eff[j] >= cut.off){
          ar.prob[skel.best] <- ar.prob[skel.best]+postprob.eff[j]
        } else {
          ar.prob[skel.best] <- ar.prob[skel.best]+0
        }
      }
      ar.prob <- ar.prob/sum(ar.prob)
      if(length(aset)==1){
        comb.best=aset
      } else {
        ifelse(sum(n)<=(cohortsize*ncohort-cohortsize),comb.best<-sample(1:ncomb,1,prob=ar.prob),comb.best<-which.max(peff.hat.aset))
        #unless it is the last sample(cohort), we keep randomizing
      }
    } #END of Strategy 4
    #################################################################################################################

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
bpocrm.sim<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb, ar.strategy){
  ncomb=length(p0)
  
  comb.select<-y<-z<-n<-matrix(nrow=ntrial,ncol=ncomb)
  nstopf=nstops=0
  
  for(i in 1:ntrial){
    result<-bpocrm(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb,ar.strategy)
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
  sim.output=list(True.tox=p0, True.eff=q0, Selection.pct=colMeans(comb.select)*100, Treat.avg=colMeans(n))
  return(sim.output)
}
##########'bpocrm.sim' end here






