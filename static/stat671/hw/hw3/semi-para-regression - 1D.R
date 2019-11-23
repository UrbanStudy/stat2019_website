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
## import the data in 1d 

library(readr)
hmw3_data1 <- read.csv("hmw3-data1.csv")

n=10
n.plot=100
p=1
x = hmw3_data1$x
y = hmw3_data1$y
x.plot=seq(from=min(x),to=max(x),length.out=n.plot)
y.plot=seq(from=min(y),to=max(y),length.out=n.plot)
plot(x.plot,y.plot,col=2,type='l')
points(x,y,pch=16)

## compute the classifier 
lambda=0.01
theta=1
I=diag(rep(1,n))
kk=outer(1:n,1:n,Vectorize(function(i,j) k(x[i],x[j])))
alpha = solve(kk + lambda*diag(rep(1,n)))%*%(y-x*theta)

## evaluate the classifier 

k.x=outer(1:n,1:n.plot,Vectorize(function(i,j) k(x[i],x.plot[j])))
hat.y=t(k.x)%*%alpha +x.plot*theta

plot(x.plot,y.plot,col=2,type='l')
lines(x.plot,hat.y)
legend("bottomleft",legend=c("true", "estimated"),
       col=c(2,1),lty=c(1,1))  




library(KSPM)
fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ y, kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))

hat.y=fitted(fit)
predict(fit, interval = "confidence")
hat.y <- predict(fit, newdata.linear = x.plot, interval = "none", level = 0.95)


equation1 <- matrix(c((kk + lambda*I),t(x)%*%kk,x,t(x)%*%x),nrow=2,ncol=2)
equation2 <- matrix(c(1,t(x)%*%x),nrow=2,ncol=1)
solve(equation1,equation2)