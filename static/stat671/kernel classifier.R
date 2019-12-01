
# 1. Simple classifier -------------------------------------------------------

# 1.1 code the simple kernelized classifier presented in class
data <- data.frame(matrix(runif(200,min = -10, max=10),nrow = 100, ncol = 2))
data$y <- ifelse(data$X1<=data$X2,1,-1)

library(ggplot2)
ggplot(data = data, aes(x = X1, y = X2, color = as.factor(y))) + 
  geom_point(size = 2) +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  theme(legend.position = "none")

kernel <- function(a,b){
  return(a %*% t(b))
}

kernel_class <- function(new,data){
  data_plus <- data[data$y==1,]
  data_minus <- data[data$y==-1,]
  
  for(i in 1:nrow(data_plus)){
    a_plus = 1/nrow(data_plus)*sum(kernel(new,data_plus[i,c(-3)]))
  }
  for(i in 1:nrow(data_minus)){
    a_minus = 1/nrow(data_minus)*sum(kernel(new,data_minus[i,c(-3)]))
  }
  
  for(i in 1:(nrow(data_plus)-1)){
    for (j in (i+1):nrow(data_plus)){
      b_plus = -0.5*(1/nrow(data_plus)^2)*sum(kernel(as.matrix(data_plus[i,c(-3)]),data_plus[j,c(-3)]))
    }
  }
  for(i in 1:(nrow(data_minus)-1)){
    for (j in (i+1):nrow(data_minus)){
      b_minus = 0.5*(1/nrow(data_minus)^2)*sum(kernel(as.matrix(data_minus[i,c(-3)]),data_minus[j,c(-3)]))
    }
  }
  
  class = a_plus+a_minus+b_plus+b_minus
  return(ifelse(class>=0,1,-1))
}


# 1.2 Show 3 examples of results obtain with simulated data in 2 dimensions, using 3 different kernels

kernel <- function(a,b){
  return(a %*% t(b))
}

kernel_class(c(1,3),data)
kernel_class(c(1,-1),data)
kernel_class(c(0,1),data)

kernel <- function(a,b){
  return((a %*% t(b) +1)^2)
}

kernel_class(c(1,3),data)
kernel_class(c(1,-1),data)
kernel_class(c(0,1),data)

kernel <- function(a,b){
  return(exp((-0.5)*(a-b)^2))
}

kernel_class(c(1,3),data)
kernel_class(c(1,-1),data)
kernel_class(c(0,1),data)

# Run this algorithm on the iris data set. Create a classifier for the lable "I.setosa" versus "I.versicolor"
iris <- read.csv("~/PDX Google Drive/LAPTOP Backup/STAT 610 Machine Learning I-III/III/iris.csv")

