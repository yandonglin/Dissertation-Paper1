rm(list=ls())

post.tox<-function(a,p,y,n){
  s2=1.34
  lik=1
  for(j in 1:length(p)){
    pj=p[j]**exp(a)
    lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
  }
  return(lik*exp(-0.5*a*a/s2));
}

#the posterior mean of ptox
posttoxf <- function(a, p, y, n, j) { p[j]^(exp(a))*post.tox(a, p, y, n); }

post.eff<-function(b,q,z,n){
  s2=1.34
  lik=1
  for(j in 1:length(q)){
    qj=q[j]**exp(b)
    lik=lik*qj^z[j]*(1-qj)^(n[j]-z[j]);
  }
  return(lik*exp(-0.5*b*b/s2));
}

#the posterior mean of peff
postefff <- function(b, q, z, n, j) {q[j]^(exp(b))*post.eff(b, q, z, n); }


### run a trial 	
ncomb = 4;   #number of combos
nord.tox=1
nord.eff=2*ncomb-1
p.skel= c(0.1,0.15,0.2,0.3)
p.skel=t(as.matrix(p.skel))
q.skel<-matrix(nrow=nord.eff,ncol=ncomb)
q.skel[1,]<-c(0.4,0.5,0.6,0.7)
q.skel[2,]<-c(0.5,0.6,0.7,0.6)
q.skel[3,]<-c(0.6,0.7,0.6,0.5)
q.skel[4,]<-c(0.7,0.6,0.5,0.4)
q.skel[5,]<-c(0.5,0.6,0.7,0.7)
q.skel[6,]<-c(0.6,0.7,0.7,0.7) 
q.skel[7,]<-c(0.7,0.7,0.7,0.7)

y=c(0,0,0,1);  #number of toxicity/responses at each dose level
z=c(0,0,1,1);   #number of efficacy at each dose level
n=c(1,1,2,2);  #number of treated patients at each dose level

ptox.hat = numeric(ncomb); # estimate of toxicity prob
peff.hat = numeric(ncomb); # estimate of efficacy prob
comb.select=rep(0,ncomb); # a vector of indicators for dose selection

mtox.sel=1 #suppose we only have one tox model, 


# efficacy model selection criteria
mprior.eff=1/nord.eff
marginal.eff = rep(0, nord.eff);
for(k in 1:nord.eff)
{
  marginal.eff[k] = integrate(post.eff,lower=-Inf,upper=Inf, q=q.skel[k,], z=z, n=n)$value;
}

################################################################################
postprob.eff = (marginal.eff*mprior.eff)/sum(marginal.eff*mprior.eff);
################################################################################

aset=c(1,2,3,4) # suppose every dose is acceptable
#round postprob.eff
#postprob.eff=c(0.229,0.172,0.089,0.065,0.196,0.137,0.112)  
###########################################################################################################
#####################Original and Alternative Randomization Strategy#######################################
###########################################################################################################
#Original: 
meff.sel=1 #suppose we select model 1
for(j in 1:ncomb){
  peff.hat[j] = integrate(postefff,lower=-Inf,upper=Inf, q.skel[meff.sel,], z, n,j)$value/marginal.eff[meff.sel]; 
}
peff.hat.aset=rep(0,ncomb)
peff.hat.aset[aset]=peff.hat[aset]

#Original Method:
  ar.prob=peff.hat.aset/sum(peff.hat.aset)
ar.prob0 <- ar.prob
 #end of original strategy


# ALternative Randomization Strategy 3:
#suppose max sample size is 30, we observed 6 so far
ncohort=30;cohortsize=1

  ar.prob <- rep(0,ncomb)
  #nord.eff is the number of candidate efficacy models
  nord.eff.limit<-ceiling(nord.eff*( ((ncohort*cohortsize-sum(n))/ncohort*cohortsize))**2  )
  #nord.eff.limit is the number of models we will consider
  if (nord.eff.limit ==0) {
    cut.off=max(postprob.eff)
  } else {
    cut.off.order <- order(postprob.eff,decreasing=F)[nord.eff-nord.eff.limit+1]
    cut.off <- postprob.eff[cut.off.order]
  }
  
  for (j in 1:length(postprob.eff)){
    skel.best <- which.max(q.skel[j,]) #The best dose by each model

    skel.best <- ifelse(is.element(skel.best,aset),skel.best,max(aset))
    
    if (postprob.eff[j] >= cut.off){
      ar.prob[skel.best] <- ar.prob[skel.best]+postprob.eff[j]
    } else {
      ar.prob[skel.best] <- ar.prob[skel.best]+0
    }
  }
  ar.prob2.raw <- ar.prob
  ar.prob <- ar.prob/sum(ar.prob)
ar.prob2 <- ar.prob
 #END of Strategy 3


# ALternative Randomization Strategy 3:
#suppose max sample size is 30, we observed 6 so far
ncohort=30;cohortsize=1

ar.prob <- rep(0,ncomb)
#nord.eff is the number of candidate efficacy models
nord.eff.limit<-ceiling(nord.eff*( ((ncohort*cohortsize-sum(n))/ncohort*cohortsize))**3  )
#nord.eff.limit is the number of models we will consider
if (nord.eff.limit ==0) {
  cut.off=max(postprob.eff)
} else {
  cut.off.order <- order(postprob.eff,decreasing=F)[nord.eff-nord.eff.limit+1]
  cut.off <- postprob.eff[cut.off.order]
}

for (j in 1:length(postprob.eff)){
  skel.best <- which.max(q.skel[j,]) #The best dose by each model
  
  skel.best <- ifelse(is.element(skel.best,aset),skel.best,max(aset))
  
  if (postprob.eff[j] >= cut.off){
    ar.prob[skel.best] <- ar.prob[skel.best]+postprob.eff[j]
  } else {
    ar.prob[skel.best] <- ar.prob[skel.best]+0
  }
}
ar.prob3.raw <- ar.prob
ar.prob <- ar.prob/sum(ar.prob)
ar.prob3 <- ar.prob
#END of Strategy 3
###############################################################################################################
ar.prob0
ar.prob2.raw
ar.prob2
ar.prob3
ar.prob3.raw







################################################################################################################
#Model selection criteria
nord.eff <- 11;ncohort <- 60;cohortsize <- 1;
n <- seq(1:59)
drop.rate <-0.5
L <- ceiling(nord.eff*(( (ncohort*cohortsize-n)/ncohort*cohortsize))**drop.rate)
windows()
par(mfrow=c(1,3))
plot(n,L,main="Drop Rate=0.5",ylab="L'", xlab="Observed Sample Size n", cex.lab=1.5)

drop.rate <-1
L <- ceiling(nord.eff*(( (ncohort*cohortsize-n)/ncohort*cohortsize))**drop.rate)
plot(n,L,main="Drop Rate=1",ylab="L'", xlab="Observed Sample Size n", cex.lab=1.5)

drop.rate <-2
L <- ceiling(nord.eff*(( (ncohort*cohortsize-n)/ncohort*cohortsize))**drop.rate)
plot(n,L,main="Drop Rate=2",ylab="L'", xlab="Observed Sample Size n", cex.lab=1.5)

#################################################################################################################


