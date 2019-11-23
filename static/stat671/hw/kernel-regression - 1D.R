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
k = function(x,y)
  return(k5(x,y))
#################################
## generate some data in 1d #####
#################################
n=10
n.plot=100
p=1
x = seq(from=0,to=2*pi,length.out=n)
x.plot=seq(from=0,to=2*pi,length.out=n.plot)
f=sin(p*x)
f.plot=sin(p*x.plot)
y = f+rnorm(n,sd=0.2)
plot(x.plot,f.plot,col=2,type='l')
points(x,y,pch=16)

##################################
## compute the classifier ########
##################################
lambda=0.01
kk=outer(1:n,1:n,Vectorize(function(i,j) k(x[i],x[j])))
alpha = solve(kk + lambda*n*diag(rep(1,n)))%*%y


##################################
## evaluate the classifier #######
## over a grid             #######
##################################
k.x=outer(1:n,1:n.plot,Vectorize(function(i,j) k(x[i],x.plot[j])))
lines(x.plot,t(k.x)%*%alpha)
legend("bottomleft",legend=c("true", "estimated"),
       col=c(2,1),lty=c(1,1))  

