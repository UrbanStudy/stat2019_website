# my first kernel classifier
rm(list=ls())
set.seed(0.1)
#library(expm)
#library(nlme)
##############################
## A few kernel functions
##############################
k1 = function(x,y)
  return(sum(x*y))

k2 = function(x,y)
  return(sum(x*y)+20)
k3 = function(x,y)
  return((1+sum(x*y))^2)
d=4
k4 = function(x,y)
  return((1+sum(x*y))^d)
sigma=10
k5 = function(x,y)
  return(exp(-sum((x-y)^2)/(2*sigma^2)))
kappa=1
theta=1
k6 = function(x,y)
  return(tanh(kappa*sum(x*y)+theta))
k7 = function(x,y)
  return(min(x,y))
k = function(x,y)
  return(k4(x,y))

##################################
## compute the kernel matrix ########
##################################
n=100
x = c(1:n)
lambda=1e-1
K=outer(1:n,1:n,Vectorize(function(i,j) k(x[i],x[j])))
E=eigen((K+lambda*diag(1,n)),symmetric=T)
K.sqrt=E$vectors%*%diag(sqrt(abs(E$values)))%*%t(E$vectors)
pdf("eigen-fourth-power-plus-offset.pdf")
plot(K.sqrt[,1],ylim=range(K.sqrt,na.rm=T),type='l',col=1)
for (i in (2:n))
  lines(K.sqrt[,i],col=i)
dev.off()
