#Getting the data frame ready
#Tested scenarios
rm(list=ls())
library(ggplot2)
#The following generates df1, which is the the dataframe for does-reposonse curves.
p <- matrix(ncol=6,nrow=4)
p[1,] <- c(0.05,0.10,0.20,0.28,0.50,0.5)  #T1 MTD=7 MODERATE T1
p[2,] <- c(0.05,0.10,0.20,0.28,0.40,0.55) #T2 MTD=4 T2
p[3,] <- c(0.05,0.10,0.15,0.20,0.35,0.40) #T3 
p[4,] <- c(0.05,0.05,0.05,0.05,0.05,0.05) #T4 MTD=6 low T4
q <- matrix(ncol=6,nrow=4)
q[1,] <- c(0.05,0.13,0.25,0.38,0.50,0.63) #R1 increasing
q[2,] <- c(0.05,0.23,0.47,0.70,0.70,0.70) #R2 plateau1
q[3,] <- c(0.7,0.7,0.7,0.7,0.7,0.7)       #R3 Flat
q[4,] <- c(0.05,0.30,0.55,0.43,0.30,0.20) #R4 Peak
Tlabel <- c('T1','T2','T3','T4')
Rlabel <- c('R1','R2','R3','R4')
dose <- Tox <- Eff <- df1 <- NULL
for (i in 1:4){
  for (j in 1:4){

    for (x in 1:6) {
      
      Tox[x] <- p[i,x]
      Eff[x] <- q[j,x]
      dose[x] <- x
      
    }
    df.temp <- data.frame(dose,Tox,Eff,Efficacy=Rlabel[j],T_Curve=Tlabel[i])
    df1 <- rbind(df1,df.temp)
  }
}
df1
df.scenario1 <- df1[,c(1,3,4,5)] #only keep efficacy rate
df.scenario1 <- cbind(df.scenario1,type="Efficacy")
df.scenario2 <- df1[,c(1,2,4,5)] #only keep toxicity rate
df.scenario2 <- cbind(df.scenario2,type="Toxicity")
colnames(df.scenario1) <- colnames(df.scenario2)<- c("Dose","Rate","Efficacy_Curve","Toxicity_Level",  "Type")
df.scenario <- rbind(df.scenario1,df.scenario2)


#Making the plot
windows()
p <- ggplot(data=df.scenario, aes(x=Dose,y=Rate,group=Type,colour=Type,lty=Type,order=Type))+geom_point() 
p <- p+scale_colour_manual(values = c("blue","red")) 
p <- p+geom_line()
p <- p+facet_grid(Toxicity_Level~Efficacy_Curve)
p <- p+labs(title="Simulated Scenarios")
p <- p+ scale_x_continuous(breaks = c(1,2,3,4,5,6))
p <- p+ scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=c('0.1','0.2','MTD 0.3','0.4','0.5','0.6','0.7'))
p <- p+geom_hline( aes(yintercept = 0.3),lty=2)
p <- p+xlab("Dose Level")+ylab("Toxicity/Efficacy Rate")
# p <- p+geom_text(aes(x=1.2,y=0.3,label='Target DLT Rate =0.3',size=0.5))
p
