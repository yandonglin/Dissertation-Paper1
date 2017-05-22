#Plotting the results
rm(list=ls())
library(xlsx)
library(xlsxjars)
library(rJava)
library(ggplot2)

#reading result data
df1 <- read.xlsx(file="C:/Users/dongl/Desktop/Dissertation Research/paper-2/confirmed result.xlsx",sheetIndex=1)
colnames(df1)[c(1,2)] <- c("Toxicity_Level","Efficacy_Level")

df1$Method[df1$Method == "1"] <- "Original"
df1$Method[df1$Method == "2"] <- "Modified"
df1$Toxicity_Level[df1$Toxicity_Level== "1"] <- "Low Toxicity"
df1$Toxicity_Level[df1$Toxicity_Level== "2"] <- "Moderate Toxicity"
df1$Efficacy_Level[df1$Efficacy_Level== "1"] <- "Increasing Efficacy"
df1$Efficacy_Level[df1$Efficacy_Level== "2"] <- "Plateau Efficacy"
df1$Efficacy_Level[df1$Efficacy_Level== "3"] <- "Peak Efficacy"
df1 <- subset(df1,Toxicity_Level!="3")

df1$Method <- factor(df1$Method,levels=c("Original","Modified"))


windows()
p <- ggplot(df1,aes(size,fill=Method))+
  geom_bar(aes(weight=selection),width=5,position = position_dodge(width=5))+facet_grid(Toxicity_Level~Efficacy_Level)
p <- p+labs(title = "Selection of Optimal Dose")+ylab("Percent")+scale_x_continuous(breaks=c(0,12,24,36,48))+xlab("Adaptive Randomization Size")

#p <- p+scale_fill_manual(values = c("grey","black"))
p



windows()
p <- ggplot(df1,aes(size,fill=Method))+
  geom_bar(aes(weight=treated),width=5, position = position_dodge(width=5))+facet_grid(Toxicity_Level~Efficacy_Level)
p <- p+labs(title = "Patient treated at Optimal Dose")+ylab("Percent")+scale_x_continuous(breaks=c(0,12,24,36,48))+xlab("Adaptive Randomization Size")
p

# add on plot: toxicity comparison
df2 <- read.xlsx(file="C:/Users/dongl/Desktop/Dissertation Research/paper-2/Confirmed result- toxicity.xlsx",sheetIndex=1)

colnames(df2)[c(1,2)] <- c("Toxicity_Level","Efficacy_Level")

df2$Method[df2$Method == "1"] <- "Original"
df2$Method[df2$Method == "2"] <- "Modified"
df2$Toxicity_Level[df2$Toxicity_Level== "1"] <- "Low Toxicity"
df2$Toxicity_Level[df2$Toxicity_Level== "2"] <- "Moderate Toxicity"
df2$Efficacy_Level[df2$Efficacy_Level== "1"] <- "Increasing Efficacy"
df2$Efficacy_Level[df2$Efficacy_Level== "2"] <- "Plateau Efficacy"
df2$Efficacy_Level[df2$Efficacy_Level== "3"] <- "Peak Efficacy"
df2 <- subset(df2,Toxicity_Level!="3")
df2$Method <- factor(df2$Method,levels=c("Original","Modified"))

windows()
p <- ggplot(df2,aes(x=size,y=Total.Tox,colour=Method))
p <- p+geom_line()
p <- p+geom_point()
p <- p+facet_grid(Toxicity_Level~Efficacy_Level)
p


windows()
boxplot(Total.Tox~Method, data=df2,main="Comparing Toxicity Response by Design Method",
        ylab="Percent of patients experiencing toxicity")

