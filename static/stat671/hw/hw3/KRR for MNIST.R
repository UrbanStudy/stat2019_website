# Load the MNIST digit recognition dataset into R
# http://yann.lecun.com/exdb/mnist/
# assume you have all 4 files and gunzip'd them
# creates train$n, train$x, train$y  and test$n, test$x, test$y
# e.g. train$x is a 60000 x 784 matrix, each row is one digit (28x28)
# call:  show_digit(train$x[5,])   to see a digit.
# brendan o'connor - gist.github.com/39760 - anyall.org
rm(list=ls())
library(pROC)
library(Matrix)
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
  train <<- load_image_file('train-images.idx3-ubyte')
  test <<- load_image_file('t10k-images.idx3-ubyte')
  
  train$y <<- load_label_file('train-labels.idx1-ubyte')
  test$y <<- load_label_file('t10k-labels.idx1-ubyte')  
}


show_digit <- function(arr784, col=gray(128:1/128)) {
  image(matrix(arr784, nrow=28)[,28:1], col=col,xlab='',ylab='')
}

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

k5 = function(x,y)
  return(exp(-sum((x-y)^2)/(2*sigma2)))        
kappa=1
theta=1
k6 = function(x,y)
  return(tanh(kappa*sum(x*y)+theta))
k = function(x,y)
  return(k5(x,y))
###############################
## load the MNIST data for once
###############################
print("loading MNIST")
load_mnist()
###############################
# reduce to digit 1 against digit 7
###############################
digit1=4
digit2=8
ind.train = which(train$y==digit1 | train$y==digit2)
ind.test= which(test$y==digit1 | test$y==digit2)
train$x=train$x[ind.train,]
train$y=ifelse(train$y[ind.train]==digit1,1,-1)
test$x=test$x[ind.test,]
test$y=ifelse(test$y[ind.test]==digit1,1,-1)
print(sprintf("digits %d and %d n= %d m = %d",digit1,digit2,dim(train$x)[1],dim(test$x)[1]))
##################################
# add a column of 1
##################################
#train$x=cbind(train$x,as.vector(rep(1,times=dim(train$x)[1])))
#test$x=cbind(test$x,as.vector(rep(1,times=dim(test$x)[1])))
#########################################################
# center the training set and apply it to the test set
#########################################################
train$x = scale(train$x,scale=F)
test$x = scale(test$x,center=attributes(train$x)$scaled,scale=F)
###################################
## sample a trainig set 
###################################
#set.seed(0.1) no set seed
n=800 #even
print(sprintf("After samplng n = %d",n))
ind1=sample(which(train$y==1),size=n/2)
ind2=sample(which(train$y==-1),size=n/2)
ind=c(ind1,ind2)
train$x=train$x[ind,]
train$y=train$y[ind]
###################################
# kernel ridge regression 
###################################
lambda = 0.001
sigma2=((range(train$x)[2]-range(train$x)[1])/5)^2
print("Computing kk")
kk=outer(1:n,1:n,Vectorize(function(i,j) k(train$x[i,],train$x[j,])))
print("computing alpha")
alpha = solve(kk + lambda*n*diag(rep(1,n)))%*%train$y
###################################
## evaluate
###################################
m=dim(test$x)[1]
print("computing kkt")
kkt=outer(1:m,1:n,Vectorize(function(i,j) k(test$x[i,],train$x[j,])))
y.hat=sign(kkt%*%alpha)
######################################
## show the results
######################################
nb.to.show=25 # a square
A=table(test$y,y.hat)
print(A)
print(sprintf("%d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat)
if(length(ind)<=nb.to.show) ind.sample=ind else ind.sample=sample(ind,size=nb.to.show)
par(mfcol=c(sqrt(nb.to.show),sqrt(nb.to.show)))
for (i in ind.sample)
  show_digit(test$x[i,])