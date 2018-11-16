# BioComputing Exercise 10

rm(list=ls())
setwd("~/Desktop/biocomp2018/exercise10/Biocomp-Fall2018-181109-Exercise10")
data=read.csv("data.txt")
library(reshape2)
library(ggplot2)

# Function for the quadratic model
quadratic = function(p,x,y){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma=exp(p[4])
  
  expected=B0+B1*x+(B2*x*x)
  nll= -sum(dnorm(x=y,mean=expected, sd=sigma, log=TRUE))
  #calculates nll for any
  return(nll)
}

qInitial=c(1,1,1,1)

qFit=optim(par = qInitial, fn=quadratic, x=data$x, y=data$y)
print(qFit)

# Function for the linear model
linear= function(p,x,y){
  C0=p[1]
  C1=p[2]
  sigma=exp(p[3])
  
  expected=C0+C1*x
  nll=-sum(dnorm(x=y, mean = expected, sd=sigma, log = TRUE))
  return(nll)
}
lInitial=c(1,1,1)
lFit=optim(par=lInitial,fn=linear, x=data$x, y=data$y)
print(lFit)

# Evaluation
testStat=2*(lFit$value-qFit$value) 
1-pchisq(q=testStat,df=1)

# Part 2

speciesModel = function(t, y, p){
  # arguments = time, state variable y, parameters p
  N1=y[1]
  N2=y[2]
  
  R1=p[1]
  a11=p[2]
  a12=p[3]
  R2=p[4]
  a22=p[5]
  a21=p[6]
  
  dN1dt = R1*(1-N1*a11-N2*a12)*N1
  
  dN2dt = R2*(1-N2*a22-N1*a21)*N2
  
  return(list(c(dN1dt,dN2dt)))
}
# a12<a11
# a21<a22

# Case 1
y0=c(2,2)
Times=1:100
params1=c(0.5,0.0095,0.002,0.5,0.0085,0.0001)
#params=(R1,a11,a12,R2,a22,a21)
sim1=ode(y=y0,times=Times,func=speciesModel,parms=params1)
out1=data.frame(time=sim1[,1],pop1=sim1[,2],pop2=sim1[,3])
out1=melt(out1,id.vars = "time")
ggplot(data=out1, aes(x=time,y=value))+geom_line(aes(color=variable))
#In this case, a11>a12 and a22>a21, and the populations coexist

# Case 2
y0=c(2,2)
Times=1:1000
params2=c(0.5,0.009,0.008,0.5,0.008,0.0095)
#params=(R1,a11,a12,R2,a22,a21)
sim2=ode(y=y0,times=Times,func=speciesModel,parms=params2)
out2=data.frame(time=sim2[,1],pop1=sim2[,2],pop2=sim2[,3])
out2=melt(out2,id.vars = "time")
ggplot(data=out2, aes(x=time,y=value))+geom_line(aes(color=variable))
# In this case, a21>a22 and population 2 goes to extinction

# Case 3
y0=c(2,2)
Times=1:100
params3=c(0.5,0.007,0.009,0.5,0.007,0.002)
#params=(R1,a11,a12,R2,a22,a21)
sim3=ode(y=y0,times=Times,func=speciesModel,parms=params3)
out3=data.frame(time=sim3[,1],pop1=sim3[,2],pop2=sim3[,3])
out3=melt(out3,id.vars = "time")
ggplot(data=out3, aes(x=time,y=value))+geom_line(aes(color=variable))
# In this case, a12>a11 and population 1 geos to extinction



