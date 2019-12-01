#Function to generate draws from the joint density of beta and sigma.sq
lm.gprior <- function(y,X,g=dim(X)[1],nu0=1,
                      s20=try(summary(lm(y~-1+X))$sigma^2,silent=TRUE),
                      S=1000){
  
  n<-dim(X)[1] ; p<-dim(X)[2]
  Hg<- (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)
  SSRg<- t(y)%*%( diag(1,nrow=n)  - Hg ) %*%y
  
  s2<-1/rgamma(S, (nu0+n)/2, (nu0*s20+SSRg)/2 )
  
  Vb<- g*solve(t(X)%*%X)/(g+1)
  Eb<- Vb%*%t(X)%*%y
  
  E<-matrix(rnorm(S*p,0,sqrt(s2)),S,p)
  beta<-t(  t(E%*%chol(Vb)) +c(Eb))
  
  list(beta=beta,s2=s2)                                
}


#load the oxygen intake data: we are interested in comparing the difference in response
# to two different treatments also accounting for age
data.oxygen<-dget("yX.o2uptake")

#run regression
regOxy.intake <- lm.gprior(y=data.oxygen[,"uptake"],
                         X=data.oxygen[,-1],
                         S=10000)

par(mfrow=c(2,2))
for(j in 1:4){
  plot(density(regOxy.intake$beta[,j]),
       xlab=bquote(beta[.(j-1)]),
       main=bquote(paste("p(",beta[.(j-1)],"|y,X,",sigma^2,")")))
  abline(v=0,lty=3,lwd=2)
}

par(mfrow=c(1,1))
plot(density(regOxy.intake$s2),
     xlab=expression(sigma^2),
     main=expression(paste("p(",sigma^2,"|y,X)")))

#To conduct the comparison of interest, note that the model is:
# Yi = beta0 + beta1*trt + beta2*age + beta3*(trt*age),
#such that:
#E[Y|trt=running,age] = beta0 + beta2*age
#E[Y|trt=aerobics,age] = (beta0+beta1) + (beta2+beta3)*age
#
#So the difference in effect can be measured by:
#E[Y|trt=aerobics,age]-E[Y|trt=running,age] = beta1 + beta3*age
ages <- 20:31
diff.mat <- matrix(NA,ncol=length(ages),nrow=dim(regOxy.intake$beta)[1])
k<-1
for(a in ages){
  diff.mat[,k] <- regOxy.intake$beta[,2]+regOxy.intake$beta[,4]*a
  k <- k+1
}

colnames(diff.mat) <- ages
boxplot(diff.mat,col="cornflowerblue",pch=20,ylab="age",xlab="diff. max Oxygen intake",
        horizontal=T)
abline(v=0,lty=3,lwd=2)


#do a similar analysis with the following dataset:
diabetes.train<-dget("http://www.stat.washington.edu/~hoff/Book/Data/data/yX.diabetes.train")

#use this dataset to validate how well the model fits the data
diabetes.test<-dget("http://www.stat.washington.edu/~hoff/Book/Data/data/yX.diabetes.test")

