
######################################
#
#Replicate Wages results and Test proposaed method
#
######################################
rm(list=ls())
rm(list=ls())
source("Redesign Adaptive Randomization.R")
#source("AltFinalDose.R") #average weight on final dose selection
#####Specify the total number of combinations
d<-6

#####Specify the number of possible toxicity orderings
s<-1   

###Specify a set of toxicity skeleton values
p.skel<-c(0.01,0.08,0.15,0.22,0.29,0.36)

#####Specify the number of possible toxicity orderings
g<-11   #efficacy
psi <- 4.6
tul<-0.33 ##toxicity upper limit 
ell<-0.05 ##efficacy lower limit
cohortsize=1 ##cohort size for each inclusion
ncohort=64  ##number of cohorts
start.comb=1 ##starting combination
n.ar=48     ##size of AR phase
ntrial=1000   ##number of simulated trials 
set.seed(580)  ##random seed
###Specifiy the possible efficacy orderings of the drug combinations
q.skel<-matrix(nrow=g,ncol=d)
q.skel[1,]<-c(0.60,0.50,0.40,0.30,0.20,0.10)  
q.skel[2,]<-c(0.50,0.60,0.50,0.40,0.30,0.20)  
q.skel[3,]<-c(0.40,0.50,0.60,0.50,0.40,0.30)  
q.skel[4,]<-c(0.30,0.40,0.50,0.60,0.50,0.40)  
q.skel[5,]<-c(0.20,0.30,0.40,0.50,0.60,0.50)
q.skel[6,]<-c(0.10,0.20,0.30,0.40,0.50,0.60)  

q.skel[7,]<-c(0.20,0.30,0.40,0.50,0.60,0.60)  
q.skel[8,]<-c(0.30,0.40,0.50,0.60,0.60,0.60)  
q.skel[9,]<-c(0.40,0.50,0.60,0.60,0.60,0.60)  
q.skel[10,]<-c(0.50,0.60,0.60,0.60,0.60,0.60)  
q.skel[11,]<-c(rep(0.60,6))  

p <- matrix(ncol=d,nrow=4)
p[1,] <- c(0.05,0.10,0.20,0.28,0.50,0.5)  #T1 MTD=7 MODERATE T1
p[2,] <- c(0.05,0.10,0.20,0.28,0.40,0.55) #T2 MTD=4 T2
p[3,] <- c(0.05,0.10,0.15,0.20,0.35,0.40) #T3
p[4,] <- c(0.05,0.05,0.05,0.05,0.05,0.05) #T4 MTD=6 low T4
q <- matrix(ncol=d,nrow=4)
q[1,] <- c(0.05,0.13,0.25,0.38,0.50,0.63) #R1 increasing
q[2,] <- c(0.05,0.23,0.47,0.70,0.70,0.70) #R2 plateau1
q[3,] <- c(0.7,0.7,0.7,0.7,0.7,0.7)       #R3 Flat
q[4,] <- c(0.05,0.30,0.55,0.43,0.30,0.23) #R4 Peak

Best.lower <- c(4,4,4,6,4,4,4,4,1,1,1,1,3,3,3,3)
Best.upper <- c(4,4,4,6,4,4,4,6,4,4,4,6,3,3,3,3)
Good.lower <- c(4,4,4,4,3,3,3,3,1,1,1,1,2,2,2,2)
Good.upper <- c(4,4,4,6,4,4,4,4,4,4,4,6,4,4,4,5)
count=0
R.curve=c("R1","R2","R3","R4")
T.curve=c("T1","T2","T3","T4")
Toxicity <- Efficacy <- Select.Best <- Select.Good <- Treat.Best <- Treat.Good <- NULL
for (j in 1:4){ #R1-R4
for (i in 1:4){ #T1-T4
  count <- count+1
    cat("i= ", i, "\n");
    cat("j= ", j, "\n");  
    
    p0<-as.vector(p[i,])
    q0<-as.vector(q[j,])

    
##simulate many trials
cat("Method=Original ", "\n");

#cat("Method=Strategy 1 ", "\n");
#bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb,ar.strategy=1)
#cat("Method=Strategy 2 ", "\n");
#bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb,ar.strategy=2)
#cat("Method=Strategy 3 ", "\n");
#bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb,ar.strategy=3)
#drop.rate <- 0.5
#cat("Method=Strategy 4 ", "\n");
#bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb,ar.strategy=4)
drop.rate <- 0.5
OBJ=bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial,start.comb,ar.strategy=4)
Toxicity[count] <- T.curve[i]
Efficacy[count] <- R.curve[j]
Select.Best[count] <- sum(OBJ$Selection.pct[seq(from=Best.lower[count],to=Best.upper[count])])
Select.Good[count] <- sum(OBJ$Selection.pct[seq(from=Good.lower[count],to=Good.upper[count])])
Treat.Best[count] <- sum(OBJ$Treat.avg[seq(from=Best.lower[count],to=Best.upper[count])])
Treat.Good[count] <- sum(OBJ$Treat.avg[seq(from=Good.lower[count],to=Good.upper[count])])
  }}

Rate05=data.frame(Toxicity, Efficacy, Select.Best, Select.Good, Treat.Best, Treat.Good, Method="Rate=0.5")
library(xlsx)
write.xlsx(Rate05, "Simulation Results/Rate05.xlsx")
