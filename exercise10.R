# BioComputing Exercise 10

rm(list=ls())
setwd("~/Desktop/biocomp2018/exercise10/Biocomp-Fall2018-181109-Exercise10")
data=read.csv("data.txt")

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
2*(107.7382-128.1356)

