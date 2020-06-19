rm(list=ls())
library(pROC)
library(Matrix)
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
## load the data 
###############################
# Import data
table_habitat <- read.csv("Willamette_habitat_features.csv")

# fix some wrong values
table_habitat[42,]$Slope <- 0.000732
table_habitat[42,]$Floodplain_elevation<- 4.214
table_habitat[41,]$Slope <- 0.000586

# Add tiny value for NA
table_habitat[is.na(table_habitat)] <- 1e-8
# Reasign the index
table_habitat <-table_habitat[order(table_habitat$RKM_2008,decreasing = F),]
table_habitat$RKM_2008 <- 1:178

# Add new columns
table_habitat <- table_habitat%>%mutate(ConnectedWet_area=AllWetArea-DisconnectedWater_Area)%>% #Creat ConnectedWet_area
  mutate(perc_1_2m=perc_2m-perc_1m)%>% #Creat pure Area_2m
  mutate(Habitat_level=as.integer(ntile(table_habitat$Habitat_area, 3))) #Creat Habitat Area level

# Change to short names
original_name <- names(table_habitat) 
names(table_habitat) <- c("No","H_A","D1_A","D2_A","D1_P","D2_P","W_m_A","W_s_A","W_a_A","L_b_A","L_v_A","W_d_A","W_ia_A","W_r_A","W_L","W_A","W_m_L","S","FE","A","W_c_A","D12_P","H_A_L")
table <- rbind(original_name[1:8],names(table_habitat)[1:8],
               original_name[9:16],names(table_habitat)[9:16],
               original_name[17:24],names(table_habitat)[17:24])
# print(xtable((table)),floating=FALSE,latex.environments=NULL,booktabs=TRUE)

# Normalize the variables
 table_habitat[,1:22] <- scale(table_habitat[,1:22], center = T, scale = T)


# Remove some variables
table_habitat_16 <- table_habitat[,c(1,2,3,4,7,8,9,10,11,13,14,15,17,18,19,20,21,22)]
glimpse(table_habitat_16)

table_habitat_category <- table_habitat[,c(23,1,4,7,8,9,10,11,13,14,15,17,18,19,20)]



###############################
# Sample training and testing set
###############################
set.seed(1)
n=120 #even
print(sprintf("Sampling n = %d",n))
ind1=sample(which(table_habitat_category$H_A_L==1),size=n/3)
ind2=sample(which(table_habitat_category$H_A_L==2),size=n/3)
ind3=sample(which(table_habitat_category$H_A_L==3),size=n/3)
ind=c(ind1,ind2,ind3)
train_y <- table_habitat_category[ind,1]
train_x <- table_habitat_category[ind,2:15]
test_y <- table_habitat_category[-ind,1]
test_x <- table_habitat_category[-ind,2:15]
table(train_y);table(test_y)

###################################
# kernel ridge regression with cv
###################################
sigma2=((range(train_x)[2]-range(train_x)[1])/5)^2
print("Computing kk")
kk=outer(1:n,1:n,Vectorize(function(i,j) k(train_x[i,],train_x[j,])))
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
    alpha = solve(a=kk.cv + lambda.seq[j]*length(train.ind)*diag(rep(1,length(train.ind))),b=train_y[train.ind])
    y.hat.cv[test.ind,j]=sign(kkt.cv%*%alpha)
    #print(kkt.cv%*%alpha)
  }
}
for (j in (1:lambda.n)){
  lambda.value[j]=sum(y.hat.cv[,j]==train_y)
  print(sprintf("lambda = %f nb good = %d out of %d",lambda.seq[j],lambda.value[j],n))
}
lambda=lambda.seq[which.max(lambda.value)]
print(lambda)
print("computing alpha")
alpha.KRR = solve(kk + lambda*n*diag(rep(1,n)))%*%train_y
###################################
## transform the data
#####################################
epsilon=1e-5
T=500
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
tmp=diag(train_y)%*%z%*%beta
J.previous = mean(g(tmp))+lambda*t(beta)%*%beta
for (t in (1:T)){
  P=diag(as.vector(g.prime(tmp)))
  J.grad=z%*%P%*%train_y/n+2*lambda*beta
  beta=beta-(1/L)*J.grad
  tmp=diag(train_y,nrow=n)%*%z%*%beta
  J = mean(g(tmp))+lambda*t(beta)%*%beta
  if (t%%100==0) print(sprintf("Step %d J = %f",t,J))
  if (abs(J.previous-J)<epsilon.stop*J) break
  J.previous=J
}
alpha.KLR=kk.eigen$vectors%*%diag(1/(sqrt(kk.eigen$values)))%*%t(kk.eigen$vectors)%*%beta
#####################################
## boosting
#####################################
epsilon.stop=1e-5
print("Boosting")
g = function(x)
  return(exp(-x))
g.prime = function(x)
  return(-exp(-x))
g.prime.L=exp(5)
L=g.prime.L*sum(z*z)/n+2*lambda
beta=z%*%alpha.KRR
tmp=diag(train_y)%*%z%*%beta
J.previous = mean(g(tmp))+lambda*t(beta)%*%beta
for (t in (1:T)){
  P=diag(as.vector(g.prime(tmp)))
  J.grad=z%*%P%*%train_y/n+2*lambda*beta
  beta=beta-(1/L)*J.grad
  tmp=diag(train_y,nrow=n)%*%z%*%beta
  J = mean(g(tmp))+lambda*t(beta)%*%beta
  if (t%%100==0) print(sprintf("Step %d J = %f",t,J))
  if (abs(J.previous-J)<epsilon.stop*J) break
  J.previous=J
}
alpha.B=kk.eigen$vectors%*%diag(1/(sqrt(kk.eigen$values)))%*%t(kk.eigen$vectors)%*%beta


# ###################################
## evaluate
###################################
m=dim(test_x)[1]
print("computing kkt")
kkt=outer(1:m,1:n,Vectorize(function(i,j) k(test_x[i,],train_x[j,])))

######################################
## show the results
######################################
y.hat=round(kkt%*%alpha.KRR,0)
y.hat[which(y.hat<1)] <- 1
y.hat[which(y.hat>3)] <- 3
nb.to.show=25 # a square
A=table(test_y,y.hat)
print(A)
print(sprintf("KRR: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))


y.hat=round(kkt%*%alpha.KLR,0)
y.hat[which(y.hat<1)] <- 1
y.hat[which(y.hat>3)] <- 3
nb.to.show=25 # a square
A=table(test_y,y.hat)
print(A)
print(sprintf("KLR: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))



y.hat=kkt%*%alpha.B
y.hat[which(y.hat<1.5)] <- 1
y.hat[which(1.5<y.hat&y.hat<2.5)] <- 2
y.hat[which(y.hat>2.5)] <- 3
nb.to.show=25 # a square
A=table(test_y,y.hat)
print(A)
print(sprintf("KB: %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))


