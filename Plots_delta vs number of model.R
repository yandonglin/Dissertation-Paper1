#Making plot for redesigning adaptive randomization stage
################################################################################################################
#Delta v.s. number of models
rm(list=ls())

nord.eff <- 7;ncohort <- 30;cohortsize <- 1;
n <- seq(1:29)
drop.rate <-c(1,2,3)
L=ceiling(nord.eff*(( (ncohort*cohortsize-n)/ncohort*cohortsize))**drop.rate[1])

L <- data.frame()

L2<-   ceiling(nord.eff*(( (ncohort*cohortsize-n)/ncohort*cohortsize))**drop.rate[2])
L3<-   ceiling(nord.eff*(( (ncohort*cohortsize-n)/ncohort*cohortsize))**drop.rate[3])

Data_L2 <- cbind(L=L2,n,rate =2)
Data_L3 <- cbind(L=L3,n,rate =3)
Data_L <- as.data.frame(rbind(Data_L2, Data_L3))
Data_L[,3] <- as.factor(Data_L[,3])
names(Data_L)[3] <- 'drop.rate'
Data_L$drop.rate <- as.character(Data_L$drop.rate)
Data_L$drop.rate[Data_L$drop.rate==2] <- "Delta=2"
Data_L$drop.rate[Data_L$drop.rate==3] <- "Delta=3"

  
#making plot

library(ggplot2)
expression(delta)
legend.title <- c(expression(delta))
windows()
p <- ggplot(Data_L ) + geom_point(aes(x = n, y = L))+facet_grid(.~ drop.rate)
p <- p+labs(title='Number of Candidate Model v.s. Sample Size')
p <- p+xlab("n: Current Sample Size")
p <- p+ylab("L': Number of Candidate Models")
p <- p+labs(shape=legend.title, colour=legend.title) 
p
#End of Delta v.s. number of models

