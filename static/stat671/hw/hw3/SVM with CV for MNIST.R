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
show_digit <- function(arr784, col=gray(128:1/128),...) {
  image(x=1:28,y=1:28,matrix(arr784, nrow=28)[,28:1], col=col,xlab='',ylab='',xaxt='n',yaxt='n')
}
##############################
## logistic functions 
##############################
sigmoid = function(x)
  return(1/(1+exp(-x)))
logistic=function(x)
  return(log(1+exp(-x)))
logistic.gradient = function(x)
  return(-sigmoid(-x))
logistic.gradient2=function(x)
  return(sigmoid(x)*sigmoid(-x))
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
digit1=4
digit2=9
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
# train$x = scale(train$x,scale=F)
# test$x = scale(test$x,center=attributes(train$x)$scaled,scale=F)
###################################
## sample a trainig set 
###################################
#set.seed(0.1) no set seed
n=500 #even
print(sprintf("Sampling n = %d",n))
ind1=sample(which(train$y==1),size=n/2)
ind2=sample(which(train$y==-1),size=n/2)
ind=c(ind1,ind2)
train$x=train$x[ind,]
train$y=train$y[ind]
###################################
# kernel ridge regression with cv
###################################
#lambda = 0.001
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
alpha = solve(kk + lambda*n*diag(rep(1,n)))%*%train$y
alpha_krr=alpha
###################################
## KLR 
###################################

epsilon=1e-5
logistic.epsilon=1e-5
T=100
logistic.loss=function(alpha)
  return(mean(logistic(diag(train$y)%*%kk%*%alpha))+lambda*t(alpha)%*%kk%*%alpha)

for (t in (1:T)){
  ll=logistic.loss(alpha)
  print(sprintf("Step %d Loss = %f",t,ll))
  u=diag(train$y)%*%kk%*%alpha
  P=diag(as.vector(logistic.gradient(u)))
  W=diag(as.vector(logistic.gradient2(u)))
  V=kk%*%W%*%kk+2*lambda*n*kk
  alpha=solve(a=V+diag(epsilon,nrow=n),b=-(kk%*%P%*%train$y+2*lambda*n*kk%*%alpha)+(V+diag(epsilon,nrow=n))%*%alpha)
  if (abs(logistic.loss(alpha)-ll)/ll<logistic.epsilon) break
}
alpha_klr=alpha

##############################
## SVM
##############################
hinge.epsilon=1e-8
y<- train$y
n <- length(y)
I_n<- diag(rep(1,n))
# lambda=1e-8           # from paper

# kernel function

### sigma=(range(train$x)[2]-range(train$x)[1])/5
### sigma=255/5
# sigma=8               # from paper

# k = function(x,y) return(exp(-sum((x-y)^2)/(2*(sigma^2))))  # use kbf kernel from paper

# kk=outer(1:n,1:n,Vectorize(function(i,j) k(train$x[i,],train$x[j,])))

# kk.eigen=eigen(kk+diag(epsilon,nrow=n),symmetric=TRUE)
# kk=kk.eigen$vectors%*%diag(kk.eigen$values)%*%t(kk.eigen$vectors)

# Hinge loss functions, Quadatic loss from 4.1
hinge=function(x) return(ifelse(x>=1,0,(1-x)^2))

# Optimization functions 
Opt = function(alpha,y,kk)
  return(mean(hinge(diag(y)%*%kk%*%alpha))+lambda*t(alpha)%*%kk%*%alpha)

# PrimalSVM function
PrimalSVM = function(kk,y,I_n) return(solve(lambda*I_n+kk)%*%y)

# if (n>100) {
#  n_2 <- n/2
#  I_n2<- diag(rep(1,n_2))
#  kk_n2=outer(1:n_2,1:n_2,Vectorize(function(i,j) k(train$x[i,],train$x[j,])))
#  y_n2<- y[(0.5*n_2+1):(1.5*n_2)]
#  alpha=PrimalSVM(kk_n2,y_n2,I_n2)
#} else {
#  indices <- 1:n
#}

# Initial alpha
 alpha=PrimalSVM(kk,y,I_n) 

 
##############################
## Indices decent
##############################

S=10
for (s in (1:S)){
# Initial sv
indices0<- indices <- which(round(alpha,13)!=0 ) # non zero components of alpha

n_sv <- length(indices)
I_sv <- diag(rep(1,n_sv))
y_sv<-train$y[indices]
x_sv<-train$x[indices,]
kk_sv <- outer(1:n_sv,1:n_sv,Vectorize(function(i,j) k(x_sv[i,],x_sv[j,])))
# Calculate alpha_sv
alpha_sv<- PrimalSVM(kk_sv,y_sv,I_sv)
# Update alpha
alpha[indices]  <- alpha_sv
alpha[-indices] <- 0            # other components of alpha=0
# Update indices
indices <- which(diag(y)%*%kk%*%alpha<1)
# Until sv has not changed
if (all(length(indices)==length(indices0))) break
print(sprintf("Step %d sv = %f",s,length(indices0)))
}

##### get indices backup: alpha>quantile(alpha,.55)|alpha<quantile(alpha,.45)
##### break criteria backup: all(length(indices)==length(indices0)) && all(indices==indices0) 


##############################
## Coodinate decent
##############################
S=500
for (s in (1:S)){
old <- Opt(alpha,y,kk)
print(sprintf("Step %d Margin = %f",s,old))
# Initial sv
indices0<- indices <- which(round(alpha,13)!=0 ) # non zero components of alpha
n_sv <- length(indices)
I_sv <- diag(rep(1,n_sv))
y_sv<-train$y[indices]
x_sv<-train$x[indices,]
kk_sv <- outer(1:n_sv,1:n_sv,Vectorize(function(i,j) k(x_sv[i,],x_sv[j,])))
# Calculate alpha_sv
alpha_sv<- PrimalSVM(kk_sv,y_sv,I_sv)
# Update alpha
alpha[indices]  <- alpha_sv
alpha[-indices] <- 0            # other components of alpha=0
# Update indices
indices <- which(diag(y)%*%kk%*%alpha<1)
 
if (abs(Opt(alpha,y,kk)-old)/old<hinge.epsilon ) break
}

##############################
## Gradient decent
##############################
for (s in (1:S)){
old <- Opt(alpha,y,kk)
print(sprintf("Step %d Margin = %f",s,old))
  
  # Initial sv
indices0<- indices <- which(round(alpha,13)!=0 ) # non zero components of alpha
n_sv <- length(indices)
I_sv <- diag(rep(1,n_sv))
y_sv<-train$y[indices]
x_sv<-train$x[indices,]
kk_sv <- outer(1:n_sv,1:n_sv,Vectorize(function(i,j) k(x_sv[i,],x_sv[j,])))
# Calculate alpha_sv
alpha_sv<- PrimalSVM(kk_sv,y_sv,I_sv)
# Gradient
gradient<-lambda*kk_sv%*%alpha_sv+kk_sv%*%I_sv%*%(kk_sv%*%alpha_sv-y_sv)
# Hessian
H<-lambda*kk_sv+kk_sv%*%I_sv%*%kk_sv
# Update alpha_sv
alpha_sv<- alpha_sv-solve(H)%*%gradient  
# Update alpha
alpha[indices] <- alpha_sv
alpha[-indices] <-0

  if (abs(Opt(alpha,y,kk)-old)/old<hinge.epsilon ) break
}


alpha_svm=alpha



  ## Decomposition method by teacher
  V <- eigen(kk)$vectors
  L <- eigen(kk+diag(epsilon,nrow=n))$values
  KK_sqrt=V %*% diag(sqrt(L)) %*% solve(V)
  beta=KK_sqrt%*%train$y/(2*n*lambda) 





###################################
## evaluate for alpha and alpha0
###################################
m=dim(test$x)[1]
print("computing kkt")
kkt=outer(1:m,1:n,Vectorize(function(i,j) k(test$x[i,],train$x[j,])))
y.hat_krr=sign(kkt%*%alpha_krr)
y.hat_klr=sign(kkt%*%alpha_klr)
y.hat_svm=sign(kkt%*%alpha_svm)
######################################
## show the results
######################################
nb.to.show=25 # a square
A=table(test$y,y.hat_krr)
print(A)
print(sprintf("KRR %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))

A=table(test$y,y.hat_klr)
print(A)
print(sprintf("KLR %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat_klr)

A=table(test$y,y.hat_svm)
print(A)
print(sprintf("SVM %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))
ind=which(test$y!=y.hat_svm)


if(length(ind)<=nb.to.show) ind.sample=ind else ind.sample=sample(ind,size=nb.to.show)
par(mfcol=c(sqrt(nb.to.show),sqrt(nb.to.show)))
for (i in ind.sample)
  show_digit(test$x[i,])
