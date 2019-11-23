# Load the MNIST digit recognition dataset into R
# http://yann.lecun.com/exdb/mnist/
# assume you have all 4 files and gunzip'd them
# creates train$n, train$x, train$y  and test$n, test$x, test$y
# e.g. train$x is a 60000 x 784 matrix, each row is one digit (28x28)
# call:  show_digit(train$x[5,])   to see a digit.
# brendan o'connor - gist.github.com/39760 - anyall.org
rm(list=ls())
library(pROC)
load_mnist <- function() {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <<- load_image_file('train-images-idx3-ubyte')
  test <<- load_image_file('t10k-images-idx3-ubyte')
  
  train$y <<- load_label_file('train-labels-idx1-ubyte')
  test$y <<- load_label_file('t10k-labels-idx1-ubyte')  
}


show_digit <- function(arr784, col=gray(128:1/128), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}
load_mnist()
###############################
# reduce to digit 1 against digit 2
###############################
digit1=6
digit2=9
ind.train = which(train$y==digit1 | train$y==digit2)
ind.test= which(test$y==digit1 | test$y==digit2)
train$x=train$x[ind.train,]
train$y=train$y[ind.train]
test$x=test$x[ind.test,]
test$y=test$y[ind.test]
##################################
# add a column of 1
##################################
train$x=cbind(train$x,as.vector(rep(1,times=dim(train$x)[1])))
test$x=cbind(test$x,as.vector(rep(1,times=dim(test$x)[1])))
##################################
# center
##################################

test$x = scale(test$x,center=as.vector(colMeans(train$x)),scale=F)
train$x = scale(train$x,scale=F)
###################################
# logistic regression 
###################################
phi = function(u)
  return(1/(1+exp(-u)))
#w = as.vector(rep(0,times=784),mode="numeric")
set.seed(3)
lambda.1=100
lambda.0=0.01
w=rnorm(784+1)
lambda=diag(c(rep(lambda.1, times=length(w)-1),lambda.0))

#plot.roc(test$y,as.vector(test$x%*%w),col=1)
hist(as.vector(test$x[test$y==digit1,]%*%w),xlim=range(as.vector(test$x%*%w)),col=rgb(0.1,0.1,0.1,0.5),main="")
hist(as.vector(test$x[test$y==digit2,]%*%w),add=T,col=rgb(0.8,0.8,0.8,0.5),main="")
T = 20 # number of steps
for (t in (1:T)){
  print(t)
  mu = phi(train$x %*% w)
  #g = t(train$x) %*% (mu-train$y) + 2*lambda*w
  delta=diag(as.vector(mu*(1-mu)))
  #h = matrix(0.0,nrow=784,ncol=784)
  #for (i in (1:(dim(train$x)[1])))
  #  h = h + mu[i]*(1-mu[i])*(train$x[i,]%o%train$x[i,])
  #s = sample(dim(train$x)[1],size=100)
  s1 = sample(which(train$y==digit1),size=200)
  s2 = sample(which(train$y==digit2),size=200)
  s = c(s1,s2)
  g = t(train$x[s,]) %*% (mu[s]-train$y[s]) + 2*lambda%*%w
  h = t(train$x[s,])%*%(delta[s,s]%*%train$x[s,]) + 2*lambda
  d = solve(h,g)
  w = w - 0.1*d
  hist(as.vector(test$x[test$y==digit1,]%*%w),xlim=range(as.vector(test$x%*%w)),col=rgb(0.1,0.1,0.1,0.5),main="")
  hist(as.vector(test$x[test$y==digit2,]%*%w),add=T,col=rgb(0.8,0.8,0.8,0.5),main="")
  #plot.roc(test$y,as.vector(test$x%*%w),add=T,col=t+1)
  #plot.roc(train$y,as.vector(train$x%*%w),add=T,col=t+1)
}
#plot.roc(test$y,as.vector(test$x%*%w),col=1,lty=2,add=T)
#print(auc(test$y,as.vector(test$x%*%w)))
show_digit(w[1:784])
#plot.roc(test$y,as.vector(test$x%*%w),add=T)
c11=sum(test$x[test$y==digit1,]%*%w<0)
c12=sum(test$x[test$y==digit2,]%*%w<0)
c21=sum(test$x[test$y==digit1,]%*%w>=0)
c22=sum(test$x[test$y==digit2,]%*%w>=0)
print(matrix(c(c11,c21,c12,c22),nrow=2))
