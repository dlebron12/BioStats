#Name:Dayanara Lebron
#Testing for sample size given power and effects desired to be seen. 

#install.packages("easypower")
library("easypower")

#Do not need to add f or eta.sq since we are aiming for small effect (small deviation target 300bp 300+/- 40bp)
main.eff1 <- list(name = "BWB", levels = 3) 
main.eff2 <- list(name = "Time", levels = 3)
main.eff3 <- list(name = "Frag_MM", levels = 3)
int.eff1<-list(name="BWB*Time",eta.sq="med")
int.eff2<-list(name="BWB*Frag_MM",eta.sq="med")
int.eff3<-list(name="Time*Frag_MM",eta.sq="med")

n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, int1=int.eff1, int2=int.eff2, int3=int.eff3)

#Other for 2 way sample size n-cell =2
ss.2way(a=8, b=4, alpha=0.05, beta=0.8, f.A=0.8, f.B=0.6, B=90) 

## Regression
data<-read.csv("2020 02 18 11H 06M Smear Analysis Result.csv", header=T)

data<-data[0:11,]
y=data[,"Avg..Size"]
time=c(20,20,20,20,18,18,18,18,15,15,15)
Frag_MM=data[,'Enzyme']
bwb=c(0,0,1.5,1.5,0,0,1.5,1.5,0,1.5,1.5)

X_matrix=data.frame(time, Frag_MM,bwb)
jpeg("rplot.jpg", width = 350, height = 350)
boxplot(y~time, main="Insert Sizes given Reaction Times", ylab="Mean Inserts", xlab="Reaction Times")
dev.off()

jpeg("insert_FragMM.jpg", width = 350, height = 350)
boxplot(y~factor(Frag_MM), main="Insert Sizes vs Frag MM", ylab="Mean Inserts", xlab="Frag MM")
dev.off()


jpeg("insert_FragMM_Time.jpg", width = 350, height = 350)
boxplot(y~Frag_MM*time, main="Insert Sizes vs Frag MM and Time", ylab="Mean Inserts",las=2)
dev.off()



jpeg("insert_FragMM_Time_bwb.jpg", width = 350, height = 350)
boxplot(y~Frag_MM*time*bwb, main="Insert Sizes vs Frag MM, BWB amd Time", ylab="Mean Inserts",las=2)
dev.off()


##ANOVA
#Appropiate if we want to regress for in between values. Need to standarize because of differences in units
model0<-lm(y~time*Frag_MM)


#Appropiate if looking per levels of factors, R-squared of 86. 
model1<-lm(y~factor(time)*factor(Frag_MM)*factor(bwb))
summary(model1)
anova(model1)


#! Calculate effect for time and Frag_MM

#We can test differences using Tukey HSD , we clearly see them so will not use.Data 2 includes the BWB, but not time. We want to see what the effect for the BWB is.

#! Calculate effect for BWB on Enz Frag, include in size calculation.

bwb_data<-read.csv("2020\ 02\ 05\ 16H\ 34M\ Smear\ Analysis\ Result.csv",header=T)
y2<-bwb_data[3:8,"Avg..Size"]
#To create more samples for exploration/mock data we can do sampling from the distribution. 
time<-rep(20,6)
Frag_MM<-rep(1,6)
bwb<-c(rep(0.5,10),rep(1,10),rep(1.5,10),rep(2,10),rep(3,10),rep(5,10))



#Too little points, lets simulate some data from distribution mean insert size 10 samples each with 2M reads?

#sapply(1:10, function(i){mean(sample(rnorm(4*10^6,y2[1],sd=10), size=2*10^6,replace=T))})
bwb_0.5<-sample(sample(rnorm(4*10^6,y2[1],sd=5), size=2*10^6,replace=T),10)
bwb_1<-sample(sample(rnorm(4*10^6,y2[2],sd=5), size=2*10^6,replace=T),10)
bwb_1.5<-sample(sample(rnorm(4*10^6,y2[3],sd=5), size=2*10^6,replace=T),10)
bwb_2.0<-sample(sample(rnorm(4*10^6,y2[4],sd=5), size=2*10^6,replace=T),10)
bwb_3.0<-sample(sample(rnorm(4*10^6,y2[5],sd=5), size=2*10^6,replace=T),10)
bwb_5.0<-sample(sample(rnorm(4*10^6,y2[6],sd=5), size=2*10^6,replace=T),10)
y3<-c(bwb_0.5, bwb_1, bwb_1.5,bwb_2.0,bwb_3.0,bwb_5.0)
bwb<-c(rep(0.5,10),rep(1,10),rep(1.5,10),rep(2,10),rep(3,10),rep(5,10))
time<-rep(20,60)

jpeg("bwb_1x.jpg", width = 350, height = 350)
plot(y3[order(y3)], main="Values of insert sizes are linear")
dev.off()
jpeg("simulated_data_bwb_1x,jpg", width = 350, height = 350)
boxplot(y3~bwb, main="Simulated values per bwb uL")
dev.off()
y3a<-c(bwb_0.5, bwb_1, bwb_1.5,bwb_2.0,bwb_3.0)
bwb2<-c(rep(0.5,10),rep(1,10),rep(1.5,10),rep(2,10),rep(3,10))
model2a<-lm(y3a~bwb2)
summary(model2a)


#####3bug_data####
#Not tested: 1/2 20 at 1.5, 1/2 20 at 2.5, 1/2 15 at 0.

insert_size=c(239,230,243,259,234,229)
Enzyme=c(1,1,0.5,0.5,0.5,0.5)
time=c(20,20,15,15,15,15)
bwb=c(0,0,1.5,1.5,2.5,2.5)
input=(rep(250,6))
#Data #3 Effect of BWB on Enz Frag slide 
y=insert_size
X=data.frame(Enzyme, time, bwb)

model3a<-lm(y~X)
summary(model3a)

model3=lm(insert_size~factor(Enzyme)*factor(time)*bwb)
