## ----echo=F, message=F, warning=F, fig.width=12, fig.height=4, fig.align='center'----
# priors
#Assumed prior: theta ~ N(mu0,t20)
mu0<-1.9 #prior mean for theta
t20<-0.95^2 #prior variance theta
#Assumed prior: lambda ~ Gamma(nu0/2,s20*nu0/2)
nu0<-1 
s20<-.01 

#data
y<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n<-length(y) ; mean.y<-mean(y) ; var.y<-var(y)


####
G<-100 ; H<-100

mean.grid<-seq(1.505,2.00,length=G) 
prec.grid<-seq(1.75,175,length=H) 

post.grid<-matrix(nrow=G,ncol=H)

for(g in 1:G) {
  for(h in 1:H) { 
    
    post.grid[g,h]<- dnorm(mean.grid[g], mu0, sqrt(t20)) *
      dgamma(prec.grid[h], nu0/2, s20*nu0/2 ) *
      prod( dnorm(y,mean.grid[g],1/sqrt(prec.grid[h])) )
  }
}

post.grid<-post.grid/sum(post.grid)

# pdf("fig6_1.pdf",height=1.75,width=5,family="Times")
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
       xlab=expression(theta), ylab=expression(lambda) )

mean.post<- apply(post.grid,1,sum)
plot(mean.grid,mean.post,type="l",xlab=expression(theta),
     ylab=expression( paste(italic("p("),
                            theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))

prec.post<-apply(post.grid,2,sum)
plot(prec.grid,prec.post,type="l",xlab=expression(lambda),
     ylab=expression( paste(italic("p("),
                            lambda,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 

CDF.theta <- cumsum(mean.post)
CDF.prec <- cumsum(prec.post)
CI.theta<- mean.grid[c(which(CDF.theta<=0.025)[1],which(CDF.theta>=0.975)[1])]
CI.prec <- prec.grid[c(which(CDF.prec<=0.025)[1],which(CDF.prec>=0.975)[1])]


## ----echo=T, message=F, warning=F----------------------------------------
set.seed(1)
S<-2000
PHI<-matrix(nrow=S,ncol=2)
PHI[1,]<-phi<-c( mean.y, 1/var.y)
mu0<-1.9  ; t20<-0.95^2
lambda0 <- 1/t20
beta<-.01 ; alpha<-1

### Gibbs sampling
for(s in 2:S) {
  
  # generate a new theta value from its full conditional
  mu.n <-  ( mu0*lambda0 + n*mean.y*phi[2] ) / ( lambda0 + n*phi[2] )
  lambda.n <- ( lambda0 + n*phi[2] )
  phi[1] <- rnorm(1, mu.n, sqrt(1/lambda.n) )
  
  # generate a new lambda value from its full conditional
  alpha.n<- alpha+n
  beta.n<- (beta + (n-1)*var.y + n*(mean.y-phi[1])^2 )
  phi[2]<- rgamma(1, alpha.n/2, beta.n/2)
  
  PHI[s,]<-phi         
}
###


## ----echo=F, message=F, warning=F, fig.width=12, fig.height=4, fig.align='center'----
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
m1<-5
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-15
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-100
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

sseq<-1:1000
sseq2<-1000:2000#c(200,1000,4)

### Gibbs based credible intervals
CI.gibbstheta <- round(quantile(PHI[sseq2,1],probs=c(0.025,0.975)),4)
CI.gibbslambda <- round(quantile(PHI[sseq2,2],probs=c(0.025,0.975)),4)

### t-test based confidence interval
n<-length(y) ; ybar<-mean(y) ; s2<-var(y)
t.CI <- ybar+qt( c(.025,.975), n-1) *sqrt(s2/n)


## ----echo=F, message=F, warning=F, fig.width=12, fig.height=4, fig.align='center'----

par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
       xlab=expression(theta), ylab=expression(tilde(sigma)^2) ,
       xlim=range(PHI[,1]),ylim=range(PHI[,2]) )
points(PHI[sseq,1],PHI[sseq,2],pch=".",cex=1.25 )

plot(density(PHI[,1],adj=2),  xlab=expression(theta),main="",
     xlim=c(1.55,2.05),
     ylab=expression( paste(italic("p("),
                            theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))

abline(v=CI.theta,lwd=2,col="gray")
abline( v= t.CI, col="black",lwd=1)

plot(density(PHI[,2],adj=2), xlab=expression(tilde(sigma)^2),main="",
     ylab=expression( paste(italic("p("),tilde(sigma)^2,"|",italic(y[1]),
                            "...",italic(y[n]),")",sep=""))) 




## ----echo=F, message=F, warning=F, fig.width=8, fig.height=7, fig.align='center'----

par(mfrow=c(2,2),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
plot(PHI[,1],  ylab=expression(theta^{(k)}),main="",pch=20,cex=0.7,
     xlab="Gibbs iteration (k)")
plot(PHI[,2],  ylab=expression(lambda^{(k)}),main="",pch=20,cex=0.7,
     xlab="Gibbs iteration (k)")
plot(cumsum(PHI[,1])/(1:S),  ylab=expression(theta^{(k)}),main="",
     type="l",col="cornflowerblue",lwd=2,xlab="Gibbs iteration (k)")
plot(cumsum(PHI[,2])/(1:S),  ylab=expression(lambda^{(k)}),main="",
     type="l",col="cornflowerblue",lwd=2,xlab="Gibbs iteration (k)")




## ----echo=T--------------------------------------------------------------
y <- c(3.4, 2.9, NA, 1.4, 3.2, 1.8, 4.6, NA, NA, NA, 2.8, NA)
n <- length(y)

c.vals <- c(1.2, 1.7, 2.0, 1.4, 0.6)
nc <- length(c.vals) #number of censored obs

#define censored set
J.set <- which(is.na(y))

#def hyperparameters
a <- b <- 1; r <- 2

#initialize z's and theta
z <- y
z[J.set] <- c.vals
theta <- rgamma(1,shape=(a+n*r),rate=(b+sum(z)))

#cdf evaluated at c.i
Fz.nc <- pgamma(c.vals,shape=r,rate=theta)
uvec <- rep(NA,length(J.set))


#Gibbs Sampelr variables
S <- 2000
gibbs.mat <- matrix(NA,ncol=(nc+1),nrow=S)

#Gibbs Sampler
for(k in 1:S){
  
  ###sample z's
  uvec <- runif(nc, min = Fz.nc, max = rep(1,nc))
  z[J.set] <- qgamma(uvec,shape=r,rate=theta)
  
  ###sample theta
  theta <-  rgamma(1,shape=(a+n*r),rate=(b+sum(z)))
  
  ###store draws
  gibbs.mat[k,] <- c(z[J.set],theta) 
}




## ----echo=F, message=F, warning=F, fig.width=6, fig.height=8, fig.align='center'----
burnin <- 1:(S/2)

par(mfrow=c(3,2),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
plot(gibbs.mat[-burnin,2],  ylab=expression(z[8]^{(k)}),main="",pch=20,cex=0.7,
     xlab="Gibbs iteration (k)")
plot(gibbs.mat[-burnin,(nc+1)],  ylab=expression(theta^{(k)}),main="",pch=20,cex=0.7,
     xlab="Gibbs iteration (k)")
plot(cumsum(gibbs.mat[-burnin,2])/(1:(S/2)),  ylab=expression(E(Z[8])),main="",
     type="l",col="cornflowerblue",lwd=2,pch=20,cex=0.7,xlab="Gibbs iteration (k)")
plot(cumsum(gibbs.mat[-burnin,(nc+1)])/(1:(S/2)),  ylab=expression(E(theta)),main="",
     type="l",col="cornflowerblue",lwd=2,pch=20,cex=0.7,xlab="Gibbs iteration (k)")
hist(gibbs.mat[-burnin,2],  xlab=expression(paste("est. density for ",z[8])),
     main="",col="cornflowerblue")
hist(gibbs.mat[-burnin,(nc+1)],  xlab=expression(paste("est. density for ",theta)),
     main="",col="cornflowerblue")
abline(v=quantile(gibbs.mat[-burnin,(nc+1)],c(0.025,0.975)),col="red",lwd=2)


## ----echo=F, message=F, warning=F, fig.width=8, fig.height=4, fig.align='center'----
###### Intro to MCMC diagnostics
mu<-c(-3,0,3)
s2<-c(.33,.33,.33)
w<-c(.45,.1,.45)

ths<-seq(-5,5,length=100)

#### MC Sampling
set.seed(1)
S<-1000
d<-sample(1:3,S, prob=w,replace=TRUE)
th<-rnorm(S,mu[d],sqrt(s2[d]))
THD.MC<-cbind(th,d)
####

### Plot MC marginal density for theta and compare to true marginal
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
ths<-seq(-6,6,length=1000)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
       expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MC[,1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )


## ----echo=F, message=F, warning=F, fig.width=8, fig.height=4, fig.align='center'----

#### MCMC sampling
th<-0
THD.MCMC<-NULL
S<-10000
set.seed(1)
for(s in 1:S) {
  d<-sample(1:3 ,1,prob= w*dnorm(th,mu,sqrt(s2))   )
  th<-rnorm(1,mu[d],sqrt(s2[d]) )
  THD.MCMC<-rbind(THD.MCMC,c(th,d) )
}

par(mfrow=c(1,2))

Smax=1000
plot(THD.MCMC[1:Smax,1],cex=0.5,main="",ylab="theta values")
lines( mu[THD.MCMC[1:Smax,2]])
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
       expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MCMC[1:Smax,1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
#lines( mu[THD.MCMC[,2]])

###


## ----echo=F, message=F, warning=F, fig.width=8, fig.height=7, fig.align='center'----

par(mfrow=c(2,2))
Smax=c(2000,5000)
plot(THD.MCMC[1:Smax[1],1],cex=0.5,main="",ylab="theta values")
lines( mu[THD.MCMC[1:Smax[1],2]])
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
       expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MCMC[1:Smax[1],1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )

plot(THD.MCMC[1:Smax[2],1],cex=0.5,main="",ylab="theta values")
lines( mu[THD.MCMC[1:Smax[2],2]])
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
       expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MCMC[1:Smax[2],1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )

# plot(THD.MCMC[1:Smax[3],1],cex=0.5,main="",ylab="theta values")
# lines( mu[THD.MCMC[1:Smax[3],2]])
# plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
#        w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
#        w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
#        expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
# hist(THD.MCMC[1:Smax[3],1],add=TRUE,prob=TRUE,nclass=20,col="gray")
# lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
#          w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
#          w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
#lines( mu[THD.MCMC[,2]])

###

