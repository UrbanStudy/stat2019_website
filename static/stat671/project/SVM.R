library("quadprog")
library(ggplot2)

iris <- read.csv("~/Google Drive File Stream/My Drive/LAPTOP Backup/STAT 610 Machine Learning I-III/III/iris.csv")[,2:6]
train_ind <- sample(seq_len(nrow(iris)), size = 120)
train <- iris[train_ind, ]
test <- iris[-train_ind, ]

# train <- iris[c(1:40,51:90,101:140),]
# test <- iris[c(41:50,91:100,141:150),]
train$y1 <-ifelse(train$Species=="I. Setosa", 1, -1)
train$y2 <-ifelse(train[,5]=="I. virginica", 1, -1)

# set the problem data and parameters
X <- as.matrix(train[,c("Petal.length", "Petal.width")])
y1 <- as.matrix(train$y1)
y2 <- as.matrix(train$y2)
n <- dim(X)[1]

####################################################################################################################
# solve QP with quadprog and the perturbance hack
# From the documentation:
# This routine implements the dual method of Goldfarb and Idnani (1982, 1983) for solving quadratic programming
# problems of the form min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.
####################################################################################################################
# quadprog solver requires that the D matrix be symmetric positive definite.  But the SVM problem is almost always
# only non-negative definite!  As a hack, we can perturb D by a small diagonal matrix and obtain positive definite
# matrix.  Choose eps a relatively small value for the diagonal perturbance.
eps <- 5e-4

# build the system matrices
Q1 <- sapply(1:n, function(i) y1[i]*t(X)[,i])
D1 <- t(Q1)%*%Q1
d <- matrix(1, nrow=n)
b0 <- rbind( matrix(0, nrow=1, ncol=1) , matrix(0, nrow=n, ncol=1) )
A <- t(rbind(matrix(y1, nrow=1, ncol=n), diag(nrow=n)))
# call the QP solver:
sol1 <- solve.QP(D1 +eps*diag(n), d, A, b0, meq=1, factorized=FALSE)
qpsol1 <- matrix(sol1$solution, nrow=n)
plot(qpsol1)

Q2 <- sapply(1:n, function(i) y2[i]*t(X)[,i])
D2 <- t(Q2)%*%Q2
A <- t(rbind(matrix(y2, nrow=1, ncol=n), diag(nrow=n)))
# call the QP solver:
sol2 <- solve.QP(D2 +eps*diag(n), d, A, b0, meq=1, factorized=FALSE)
qpsol2 <- matrix(sol2$solution, nrow=n)
plot(qpsol2)

# predict test set
findwb <- function(a, y, X){
  nonzero <-  abs(a) > 1e-5
  W <- rowSums(sapply(which(nonzero), function(i) a[i]*y[i]*X[i,]))
  b <- mean(sapply(which(nonzero), function(i) X[i,]%*%W- y[i]))
  return(c(W,b))
}

w1 <- as.matrix(findwb(qpsol1, y1, X)[1:2])
b1 <-  findwb(qpsol1, y1, X)[3]
test$pred_setona <- as.matrix(test[,3:4])%*%w1 -b1

w2 <- as.matrix(findwb(qpsol2, y1, X)[1:2])
b2 <-  findwb(qpsol2, y1, X)[3]
test$pred_virg <- as.matrix(test[,3:4])%*%w2 -b2
test$pred <- ifelse(test$pred_setona>0, "I. Setosa",ifelse(test$pred_virg<0, "I. virginia","I. versicolor"))

# construct confusion matrix
table(test[,5], test$pred)

# plot decision bundary
findLine <- function(a, y, X){
  nonzero <-  abs(a) > 1e-5
  W <- rowSums(sapply(which(nonzero), function(i) a[i]*y[i]*X[i,]))
  b <- mean(sapply(which(nonzero), function(i) X[i,]%*%W- y[i]))
  slope <- -W[1]/W[2]
  intercept <- b/W[2]
  return(c(intercept,slope))
}

qpline1 <- findLine(qpsol1, y1, X)
qpline2 <- findLine(qpsol2, y2, X)

ggplot(train, aes(x=Petal.length, y=Petal.width)) + 
  ggtitle("Solving the SVM") +
  geom_point(aes(fill=Species), size=3, pch=21) +
  geom_abline(intercept=qpline1[1], slope=qpline1[2], size=1, aes(color="quadprog"), show.legend =TRUE) +
  geom_abline(intercept=qpline2[1], slope=qpline2[2], size=1, aes(color="quadprog"), show.legend =TRUE)



# Kernelized --------------------------------------------------------------

kernel1 <- function(a,b) {return(a%*%t(b))} # linear dot product, same result as above
kernel2 <-  function(a,b){return((a %*% t(b) +1)^6)}
k <- kernel2

D1 <- y1%*%t(y1)*k(X,X)
d <- matrix(1, nrow=n)
b0 <- rbind( matrix(0, nrow=1, ncol=1) , matrix(0, nrow=n, ncol=1) )
A <- t(rbind(matrix(y1, nrow=1, ncol=n), diag(nrow=n)))
# call the QP solver:
sol1 <- solve.QP(D1 +eps*diag(n), d, A, b0, meq=1, factorized=FALSE)
qpsol1 <- matrix(sol1$solution, nrow=n)
plot(qpsol1)

D2 <- y2%*%t(y2)*k(X,X)
A <- t(rbind(matrix(y2, nrow=1, ncol=n), diag(nrow=n)))
# call the QP solver:
sol2 <- solve.QP(D2 +eps*diag(n), d, A, b0, meq=1, factorized=FALSE)
qpsol2 <- matrix(sol2$solution, nrow=n)
plot(qpsol2)

# which(as.data.frame(qpsol2)>500)
# matrix_to_sum <- matrix(NA, nrow = 120, ncol = 2)
# for (j in 1:120) matrix_to_sum[j,] <- rowSums(as.matrix(rowSums(qpsol2[2]*y2[2]*k(as.matrix(X[j,]),as.matrix(X[2,])))))
# a1 <- as.matrix(colSums(matrix_to_sum))
# for (j in 1:120) matrix_to_sum[j,] <- rowSums(as.matrix(rowSums(qpsol2[57]*y2[57]*k(as.matrix(X[j,]),as.matrix(X[57,])))))
# a2 <- as.matrix(colSums(matrix_to_sum))
# for (j in 1:120) matrix_to_sum[j,] <- rowSums(as.matrix(rowSums(qpsol2[112]*y2[112]*k(as.matrix(X[j,]),as.matrix(X[112,])))))
# a3 <- as.matrix(colSums(matrix_to_sum))

# w <- a1 + a2 + a3
# b <- mean(y2[c(2,57,112)]- as.matrix(rbind(t(a1),t(a2),t(a3)))%*%as.matrix(w))

# test$pred_virg <- as.matrix(test[,3:4])%*%w -b
# test$pred <- ifelse(test$pred_setona>0, "I. Setosa",ifelse(test$pred_virg<0, "I. virginia","I. versicolor"))

# cm2 <- table(test[,5], test$pred)
