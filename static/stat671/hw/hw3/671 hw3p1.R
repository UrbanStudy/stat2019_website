# kernel regression
rm(list=ls())
#set.seed(0.1)
##############################
## A few kernel functions
##############################
k1 = function(x,y)
  return(sum(x*y))
k12 = function(x,y)
  return(sum(x*y)+1)
k2 = function(x,y)
  return(sum(x*y)^2)
k3 = function(x,y)
  return((1+sum(x*y))^2)
d=4
k4 = function(x,y)
  return((1+sum(x*y))^d)
sigma=1
k5 = function(x,y)
  return(exp(-sum((x-y)^2)/(2*sigma^2)))        
kappa=1
theta=1
k6 = function(x,y)
  return(tanh(kappa*sum(x*y)+theta))
kernel = function(x,y)
  return(k1(x,y))
#################################
## generate some data in 1d #####
#################################

n.p=40
n.m=40
n=n.p+n.m
library(mvtnorm)
x.p=rmvnorm(n=n.p,mean=c(1,1)+c(2,2),sigma=diag(rep(1,2)))
x.m=rmvnorm(n=n.m,mean=c(-1,-1)+c(2,2),sigma=diag(rep(2,2)))
x=rbind(x.p,x.m)
#y <- ifelse(x[,1]<=x[,2],1,-1)
y<-ifelse(x[,1]<=x[,2],1,-1)
##################################
## compute the classifier ########
##################################
lambda=0.00000001
N = nrow(x)
ident.N = diag(rep(1,N))
KK <- matrix(rep(0,N^2),N,N)
KK=outer(1:n,1:n,Vectorize(function(i,j) kernel(x[i,],x[j,])))
alpha = solve(KK + lambda*N*diag(rep(1,N)))%*%y
f=t(KK)%*%alpha
result<-ifelse(f>=0,1,-1)
y<-as.factor(y)
result<-as.factor(result)
confusionMatrix(result, y)

##################################
## evaluate the classifier #######
## over a grid             #######
##################################






