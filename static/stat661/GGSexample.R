GGSexample <- function(ksi, lambda, delta = 1e-006, n, N) 
{
  ## This is to apply a generalized Gauss-Seidel (Coordinate Ascent) algorithm to obtain the MLE of the mixing probability and the mean for Binomial/Poisson  mixture model
  ##  ##################################################################  
  ## "ksi" is the unknown mixing probability parameter  
  ## "lambda" is the unknown mean parameter for a Poisson distribution  
  ## "delta" is a tolerance limit for convergence; A default value is .000001  
  ## "n" is a column matrix, consisting of the observed frequencies of categories.   
  ## That is, the number of widows who had a specific number of children ranging   
  ## from 0 to 6 in this example  ## "N" is the total sample size  
  ## Author: Jong Sung Kim, 11/25/19  ##################################################################  
nzero <- n[1, 1] # no of widows who had no children  
size <- nrow(n) # no of categories  
index <- 1:(size - 1)  
repeat {   newksi <- (nzero * exp(lambda) - N)/(N * (exp(lambda) - 1))   
repeat {   newlambda <- (sum(index * n[2:size,  ]) * (newksi * exp(lambda) + (1 - newksi)))/((N - nzero) * newksi * exp(lambda) + N * (1 -  newksi))    
conv <- abs(newlambda - lambda) 
# absolute difference between 
# two consecutive iterates
if(conv < delta) break    
lambda <- newlambda   }
newlambda   
x <- c(newksi, newlambda)   
y <- c(ksi, lambda)   
conv <- dist(rbind(x,y)) # absolute difference
## between two consecutive iterates   
if(conv < delta) break   
ksi <- newksi   
lambda <- newlambda  }  
print("ksihat and lambdahat are:")  
return(c(newksi, newlambda)) 
} 
  
GGSexample(.75,.4,n=matrix(c(3062,587,284,103,33,4,2),7,1),N=4075) 