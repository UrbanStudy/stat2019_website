# Load the MNIST digit recognition dataset into R
# http://yann.lecun.com/exdb/mnist/
# assume you have all 4 files and gunzip'd them
# creates train$n, train$x, train$y  and test$n, test$x, test$y
# e.g. train$x is a 60000 x 784 matrix, each row is one digit (28x28)
# call:  show_digit(train$x[5,])   to see a digit.
# brendan o'connor - gist.github.com/39760 - anyall.org
rm(list=ls())
#library(pROC)
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


show_digit <- function(arr784, col=gray(128:1/128),...) {
  image(x=1:28,y=1:28,matrix(arr784, nrow=28)[,28:1], col=col,xlab='',ylab='',xaxt='n',yaxt='n')
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
  return(k3(x,y))
###############################
## load the MNIST data for once
###############################
print("loading MNIST")
load_mnist()
###############################
# reduce to digit 1 against digit 7
###############################
digit1=1
digit2=7
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
n=300 #even
print(sprintf("Sampling n = %d",n))
ind1=sample(which(train$y==1),size=n/2)
ind2=sample(which(train$y==-1),size=n/2)
ind=c(ind1,ind2)
train$x=train$x[ind,]
train$y=train$y[ind]

###################################
# kernel ridge regression with cv
###################################
sigma2=((range(train$x)[2]-range(train$x)[1])/5)^2
print("Computing kk")
kk=outer(1:n,1:n,Vectorize(function(i,j) k(train$x[i,],train$x[j,])))
lambda.seq.log10=(-5:5)
lambda.seq=10^(lambda.seq.log10)
lambda.n=length(lambda.seq)
n.fold=10
y.hat.cv=matrix(nrow=n,ncol=lambda.n)
lambda.value=vector(length=lambda.n)
which.fold=sample.int(n=n.fold,size=n,replace=T)
for (i in (1:n.fold)){
  train.ind=which(which.fold!=i)
  test.ind=which(which.fold==i)
  kk.cv=kk[train.ind,train.ind]
  kkt.cv=kk[test.ind,train.ind]
  for (j in (1:lambda.n)){
    alpha = solve(a=kk.cv + lambda.seq[j]*length(train.ind)*diag(rep(1,length(train.ind))),b=train$y[train.ind])
    y.hat.cv[test.ind,j]=sign(kkt.cv%*%alpha)
    #print(kkt.cv%*%alpha)
  }
}
for (j in (1:lambda.n)){
  lambda.value[j]=sum(y.hat.cv[,j]==train$y)
  print(sprintf("lambda = %f nb good = %d out of %d",lambda.seq[j],lambda.value[j],n))
}
lambda=lambda.seq[which.max(lambda.value)]
print(lambda)
print("computing alpha")
alpha.KRR = solve(kk + lambda*n*diag(rep(1,n)))%*%train$y
###################################
## transform the data
#####################################
epsilon=1e-5
T=5000
kk.eigen=eigen(kk+diag(epsilon,nrow=n),symmetric=TRUE)
z=kk.eigen$vectors%*%diag(sqrt(kk.eigen$values))%*%t(kk.eigen$vectors)
print("Computing the Gram matrix")
Gram=t(z)%*%z
#####################################
## logistic regression
#####################################
print("Logistic regression")
epsilon.stop=1e-5
sigmoid = function(x)
  return(1/(1+exp(-x)))
g = function(x)
  return(log(1+exp(-x)))
g.prime = function(x)
  return(-sigmoid(-x))
g.prime.L=.25
L=g.prime.L*sum(z*z)/n+2*lambda
beta=z%*%alpha.KRR
tmp=diag(train$y)%*%z%*%beta
J.previous = mean(g(tmp))+lambda*t(beta)%*%beta
for (t in (1:T)){
  P=diag(as.vector(g.prime(tmp)))
  J.grad=z%*%P%*%train$y/n+2*lambda*beta
  beta=beta-(1/L)*J.grad
  tmp=diag(train$y,nrow=n)%*%z%*%beta
  J = mean(g(tmp))+lambda*t(beta)%*%beta
  if (t%%100==0) print(sprintf("Step %d J = %f",t,J))
  if (abs(J.previous-J)<epsilon.stop*J) break
  J.previous=J
}
alpha.KLR=kk.eigen$vectors%*%diag(1/(sqrt(kk.eigen$values)))%*%t(kk.eigen$vectors)%*%beta
#####################################
## boosting
#####################################
epsilon.stop=1e-8
print("Boosting")
g = function(x)
  return(exp(-x))
g.prime = function(x)
  return(-exp(-x))
g.prime.L=exp(5)
L=g.prime.L*sum(z*z)/n+2*lambda
beta=z%*%alpha.KRR
tmp=diag(train$y)%*%z%*%beta
J.previous = mean(g(tmp))+lambda*t(beta)%*%beta
for (t in (1:T)){
  P=diag(as.vector(g.prime(tmp)))
  J.grad=z%*%P%*%train$y/n+2*lambda*beta
  beta=beta-(1/L)*J.grad
  tmp=diag(train$y,nrow=n)%*%z%*%beta
  J = mean(g(tmp))+lambda*t(beta)%*%beta
  if (t%%100==0) print(sprintf("Step %d J = %f",t,J))
  if (abs(J.previous-J)<epsilon.stop*J) break
  J.previous=J
}
alpha.B=kk.eigen$vectors%*%diag(1/(sqrt(kk.eigen$values)))%*%t(kk.eigen$vectors)%*%beta

##########################################################
###############   4.1  newton method  ####################
##########################################################
print("SVM by newton method")
# Hinge functions 
hinge=function(x)
  return(ifelse(x>=1,0,(1-x)^2))

# Hinge loss functions 
hinge.loss=function(alpha)
  return(mean(hinge(diag(train$y)%*%kk%*%alpha))+lambda*t(alpha)%*%kk%*%alpha)

## Quadatic loss by article
n <- length(train$y)
value=50
kk=outer(1:n,1:n,Vectorize(function(i,j) k(train$x[i,],train$x[j,])))
#Beta=solve(kk + lambda*n*diag(rep(1,n)))%*%train$y
Beta=z%*%alpha.KRR


while (value >=epsilon.stop*J)  {
  
  old<-(mean(hinge(diag(train$y)%*%kk%*%Beta))+lambda*t(Beta)%*%kk%*%Beta)
  print(sprintf("Step %d J = %f",t,old))
  
  
  # select support vector based on yif(xi)
  indices1 <- indices <- which(diag(train$y)%*%kk%*%Beta <1)
  indices0 <- indices <- which(diag(train$y)%*%kk%*%Beta >=1)
  
  #generate I0
  diag_I0=matrix(0,nrow=n,ncol=1)
  diag_I0<-as.vector(diag_I0)
  num_supVec_I<-length(indices1)
  num_other_I<-length(indices0)
  supVec_I<-diag_I0[indices1]<-rep(1,num_supVec_I)
  other_I<-diag_I0[indices0]<-rep(0,num_other_I)
  diag_I0=c(c(supVec_I),c(other_I))
  I0<-diag(diag_I0)
  I0<-I0+0.000001*(diag(rep(1,n)))
  
  #reorder train$y
  supVec_y<-train$y[indices1]
  other_y<-train$y[indices0]
  train$y<-c(c(supVec_y),c(other_y))
  #reorder train$x
  supVec_x<-train$x[indices1,]
  other_x<-train$x[indices0,]
  train$x<-rbind(supVec_x,other_x)
  # regenerate K
  kk=outer(1:n,1:n,Vectorize(function(i,j) k(train$x[i,],train$x[j,])))
  
  #gradient is equation(14) , a 500*1 matrix
  gradient<-2*(lambda*kk%*%Beta+kk%*%I0%*%(kk%*%Beta-train$y))
  # H = 2(lambda*k+k*I0*k)  equation (15), a 500*500 matrix
  H<-2*(lambda*kk+kk%*%I0%*%kk)
  invH<-solve(H)  
  
  # Newton method iteration
  gamma<-1
  Beta<- Beta-gamma*invH%*%gradient
  new<-(mean(hinge(diag(train$y)%*%kk%*%Beta))+lambda*t(Beta)%*%kk%*%Beta)
  value<-abs(old-new)
}

# Equation 16
Beta1<-solve((lambda*kk+kk%*%I0%*%kk))%*%kk*I0%*%diag(train$y)

# ###################################
## evaluate
###################################
m=dim(test$x)[1]
print("computing kkt")
kkt=outer(1:m,1:n,Vectorize(function(i,j) k(test$x[i,],train$x[j,])))

######################################
## show the results
######################################
y.hat=sign(kkt%*%alpha.KRR)
nb.to.show=25 # a square
A=table(test$y,y.hat)
print(A)
print(sprintf("KRR: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat)
if(length(ind)<=nb.to.show) ind.sample=ind else ind.sample=sample(ind,size=nb.to.show)
par(mfcol=c(sqrt(nb.to.show),sqrt(nb.to.show)))
for (i in ind.sample)
  show_digit(test$x[i,])

y.hat=sign(kkt%*%alpha.KLR)
nb.to.show=25 # a square
A=table(test$y,y.hat)
print(A)
print(sprintf("KLR: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat)
if(length(ind)<=nb.to.show) ind.sample=ind else ind.sample=sample(ind,size=nb.to.show)
par(mfcol=c(sqrt(nb.to.show),sqrt(nb.to.show)))
for (i in ind.sample)
  show_digit(test$x[i,])

y.hat=sign(kkt%*%alpha.B)
nb.to.show=25 # a square
A=table(test$y,y.hat)
print(A)
print(sprintf("KB: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat)
if(length(ind)<=nb.to.show) ind.sample=ind else ind.sample=sample(ind,size=nb.to.show)
par(mfcol=c(sqrt(nb.to.show),sqrt(nb.to.show)))
for (i in ind.sample)
  show_digit(test$x[i,])

y.hat=sign(kkt%*%Beta)
nb.to.show=25 # a square
A=table(test$y,y.hat)
print(A)
print(sprintf("SVM: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat)
if(length(ind)<=nb.to.show) ind.sample=ind else ind.sample=sample(ind,size=nb.to.show)
par(mfcol=c(sqrt(nb.to.show),sqrt(nb.to.show)))
for (i in ind.sample)
  show_digit(test$x[i,])

