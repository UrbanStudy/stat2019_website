rm(list=ls())
library(pROC)
library(Matrix)
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

##############################
## SVM
##############################
hinge.epsilon=1e-8
n <- length(train_y)
I_n<- diag(rep(1,n))
lambda=1e-4           # from paper
#sigma2=0.001         # from paper
sigma2=((range(train_x)[2]-range(train_x)[1])/5)^2
# kernel 
kk=outer(1:n,1:n,Vectorize(function(i,j) k(train_x[i,],train_x[j,])))

# Hinge loss functions, Quadatic loss from 4.1
hinge=function(x) return(ifelse(x>=1,0,(1-x)^2))

# Optimization functions 
Opt = function(alpha,train_y,kk)
  return(sum(hinge(diag(train_y)%*%kk%*%alpha))+lambda*t(alpha)%*%kk%*%alpha)

# PrimalSVM function
PrimalSVM = function(kk,train_y,I_n) return(solve(lambda*I_n+kk)%*%train_y)


# Initial alpha
 alpha=PrimalSVM(kk,train_y,I_n) 

 ##############################
 ## 1. Indices decent
 ##############################
 ptm <- proc.time()
 S=10
 for (s in (1:S)){
   # Initial sv
   indices0<- indices <- which(abs(alpha)<10) # non zero components of alpha
   
   n_sv <- length(indices)
   I_sv <- diag(rep(1,n_sv))
   y_sv<-train_y[indices]
   x_sv<-train_x[indices,]
   kk_sv <- outer(1:n_sv,1:n_sv,Vectorize(function(i,j) k(x_sv[i,],x_sv[j,])))
   # Calculate alpha_sv
   alpha_sv<- PrimalSVM(kk_sv,y_sv,I_sv)
   # Update alpha
   alpha[indices]  <- alpha_sv
   alpha[-indices] <- 0            # other components of alpha=0
   # Update indices
   indices <- which(diag(train_y)%*%kk%*%alpha<1)
   # Until sv has not changed
   print(sprintf("Step %d sv = %f",s,length(indices0)))
   if (all(length(indices)==length(indices0))) break   
 }
 Indices_decent <-proc.time() - ptm
 alpha_svm_indices=alpha
 
##############################
## 2. Margin Ascend
##############################
 # Initial alpha
 alpha=PrimalSVM(kk,train_y,I_n) 
 ptm <- proc.time() 
S=500
for (s in (1:S)){
old <- Opt(alpha,train_y,kk)
print(sprintf("Step %d Margin = %f",s,old))
# Initial sv
indices0<- indices <- which(round(alpha,13)!=0 ) # non zero components of alpha
n_sv <- length(indices)
I_sv <- diag(rep(1,n_sv))
y_sv<-train_y[indices]
x_sv<-train_x[indices,]
kk_sv <- outer(1:n_sv,1:n_sv,Vectorize(function(i,j) k(x_sv[i,],x_sv[j,])))
# Calculate alpha_sv
alpha_sv<- PrimalSVM(kk_sv,y_sv,I_sv)
# Update alpha
alpha[indices]  <- alpha_sv
alpha[-indices] <- 0            # other components of alpha=0
# Update indices
indices <- which(diag(train_y)%*%kk%*%alpha<1)
 
if (abs(Opt(alpha,train_y,kk)-old)/old<hinge.epsilon ) break
}
Margin_Ascend <-proc.time() - ptm
alpha_svm_margin=alpha
 
 
##############################
## Loss decent
##############################

# Initial alpha
alpha=PrimalSVM(kk,train_y,I_n) 
ptm <- proc.time() 
for (s in (1:S)){
  # Initial sv
indices0<- indices <- which(round(alpha,13)!=0  ) # non zero components of alpha
n_sv <- length(indices)
I_sv <- diag(rep(1,n_sv))
y_sv<-train$y[indices]
x_sv<-train$x[indices,]
kk_sv <- outer(1:n_sv,1:n_sv,Vectorize(function(i,j) k(x_sv[i,],x_sv[j,])))
# Calculate alpha_sv
alpha_sv<-alpha[indices]

old <- Opt(alpha_sv,y_sv,kk_sv)
print(sprintf("Step %d Loss = %f",s,old))

# Gradient
gradient<-lambda*kk_sv%*%alpha_sv+kk_sv%*%I_sv%*%(kk_sv%*%alpha_sv-y_sv)
# Hessian
H<-lambda*kk_sv+kk_sv%*%I_sv%*%kk_sv
# Update alpha_sv
alpha_sv<- alpha_sv-solve(H)%*%gradient  

  if (abs(Opt(alpha_sv,y_sv,kk_sv)-old)/old<hinge.epsilon ) break

# Update alpha
alpha[indices] <- alpha_sv
alpha[-indices] <-0
}
Loss_decent <-proc.time() - ptm

alpha_svm_loss=alpha


###################################
## evaluate for alpha and alpha0
###################################
m=dim(test_x)[1]
print("computing kkt")
kkt=outer(1:m,1:n,Vectorize(function(i,j) k(test_x[i,],train_x[j,])))

y.hat_svm_indices=kkt%*%alpha_svm_indices
y.hat_svm_indices[which(y.hat_svm_indices<1.5)] <- 1
y.hat_svm_indices[which(1.5<y.hat_svm_indices&y.hat_svm_indices<2.5)] <- 2
y.hat_svm_indices[which(y.hat_svm_indices>2.5)] <- 3

y.hat_svm_margin=kkt%*%alpha_svm_margin
y.hat_svm_margin[which(y.hat_svm_margin<1.5)] <- 1
y.hat_svm_margin[which(1.5<y.hat_svm_margin&y.hat_svm_margin<2.5)] <- 2
y.hat_svm_margin[which(y.hat_svm_margin>2.5)] <- 3

y.hat_svm_loss=kkt%*%alpha_svm_loss
y.hat_svm_loss[which(y.hat_svm_loss<1.5)] <- 1
y.hat_svm_loss[which(1.5<y.hat_svm_loss&y.hat_svm_loss<2.5)] <- 2
y.hat_svm_loss[which(y.hat_svm_loss>2.5)] <- 3

######################################
## show the results
######################################

A=table(test_y,y.hat_svm_indices)
print(A)
print(sprintf("SVM %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))

A=table(test_y,y.hat_svm_margin)
print(A)
print(sprintf("SVM %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))

A=table(test_y,y.hat_svm_loss)
print(A)
print(sprintf("SVM %d correct out of %d, or %.2f percent",sum(diag(A)),sum(A),100*sum(diag(A))/sum(A)))

time <- rbind(Indices_decent,Margin_Ascend, Loss_decent)[,3]
barplot(time, main="running time of different methods",
        ylab="second", col=c("darkgreen","GreenYellow","Goldenrod"),
        legend = rownames(time), beside=TRUE)

summary(kkt%*%alpha_svm_indices)
