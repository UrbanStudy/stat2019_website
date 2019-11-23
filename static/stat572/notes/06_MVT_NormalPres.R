## ----echo=F, message=F, eval=T, warning=F, fig.width=4.2, fig.height=2, fig.align='center'----
### sample from the multivariate normal distribution
rmvnorm<-function(n,mu,Sigma) 
{
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 )
  {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  res
}
###


### sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
     Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
     S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}
###


### reading data
Y <- dget("http://www2.stat.duke.edu/~pdh10/FCBS/Inline/Y.reading")

n <- nrow(Y); p<- ncol(Y)
y.bar <- colMeans(Y)

#initialize Sigma at sample covariance
Sigma <- cov(Y)

##prior parameters
#for theta
mu0 <- c(50,50)
L0 <- matrix(c(625,312.5,312.5,625),ncol=2)

#for Sigma
nu0 <- 4
S0 <- L0


set.seed(123)

S <- 10000
##Storage objects
theta.Mat <- matrix(NA,ncol=2,nrow=S)
Sigma.Mat <- matrix(NA,ncol=4,nrow=S)

YS<-NULL

S0.inv <- solve(S0)

for(i in 1:S){
  ##update theta
  Ln <- solve(S0.inv+n*solve(Sigma))
  mun <- Ln %*% (S0.inv %*%mu0 + n*solve(Sigma)%*%y.bar)
  theta <- c(rmvnorm(1,mun,Ln))
  
  ##update Sigma
  Sn <- S0 + (t(Y)-c(theta))%*%t(t(Y)-c(theta))
  Sigma <- solve(rwish(1,nu0+n,solve(Sn)))
  
  ### save predicted draws
  YS<-rbind(YS,rmvnorm(1,theta,Sigma))
  
  ###Save new draws
  theta.Mat[i,] <- theta
  Sigma.Mat[i,] <- c(Sigma)
}

#estimate P(theta2>theta1 | y)
pmu <- mean(theta.Mat[,2]>theta.Mat[,1])
pys <- mean(YS[,2]>YS[,1])

par(mfrow=c(1,2),mgp=c(1.75,.75,0),mar=c(3,3,1,1))
plot(x=theta.Mat[,1],y=theta.Mat[,2],xlab=expression(theta[1]),
     ylab=expression(theta[2]), cex=0.3)
abline(0,1,col="red",lwd=2)
plot(YS,xlab=expression(italic(y[1])),ylab=expression(italic(y[2])), 
     xlim=c(0,100),ylim=c(0,100), cex=0.3)
abline(0,1,col="red",lwd=2)
points(Y[,1],Y[,2],pch=16,cex=.7,col="cornflowerblue")


