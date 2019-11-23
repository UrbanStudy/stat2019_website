#### Sparrow data
load("sparrows.RData")  
fledged<-sparrows[,1] ; age<-sparrows[,2] ; age2<-age^2



#### Figure 10.1 
pdf("fig10_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(fledged~as.factor(age),range=0,xlab="age",ylab="offspring",
      col="gray")
summary(glm(fledged~age+age2,family="poisson"))
dev.off()


#### Grid-based posterior approximation
p<-3
beta0<-rep(0,p)
S0<-diag( rep(100,3))
gs<-100
LPB<-array(0,dim=rep(gs,p))

beta1<-seq(.27-1.75,.27+1.75,length=gs)
beta2<-seq(.68-1.5,.68+1.5,length=gs)
beta3<-seq(-.13-.25,-.13+.25,length=gs)

beta1<-seq(.27-2.5,.27+2.5,length=gs)
beta2<-seq(.68-2,.68+2,length=gs)
beta3<-seq(-.13-.5,-.13+.5,length=gs)



for(i in 1:gs) { for(j in 1:gs) { for(k in 1:gs) {
  theta<-beta1[i]+beta2[j]*age+beta3[k]*age^2
  LPB[i,j,k]<-dnorm(beta1[i],beta0[1],sqrt(S0[1,1]),log=TRUE)  +
              dnorm(beta2[j],beta0[2],sqrt(S0[2,2]),log=TRUE)  +
              dnorm(beta3[k],beta0[3],sqrt(S0[3,3]),log=TRUE)  +
              sum( dpois(fledged,exp(theta),log=TRUE )  )
  }} 
 cat(i,"\n") }

PB<-exp( LPB - max(LPB) )
PB<-PB/sum(PB)

PB1<-apply(PB,1,sum)
PB2<-apply(PB,2,sum)
PB3<-apply(PB,3,sum)
PB23<-apply(PB,c(2,3),sum)


## Simulation from grid approximation
S<-50000
BETAg<-matrix(nrow=S,ncol=3)
for(s in 1:S) {
i<-sample(1:gs,1,prob=PB2)
j<-sample(1:gs,1,prob=PB23[i,] )
k<-sample(1:gs,1,prob=PB[,i,j] )
BETAg[s,]<-c(beta1[k],beta2[i],beta3[j])  }


#### Figure 10.2
pdf("fig10_2.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(beta2,PB2*length(beta2)/(max(beta2)-min(beta2)) ,type="l",xlab=expression(beta[2]),ylab=expression(paste(italic("p("),beta[2],"|",italic("y)"),sep="") ) )
plot(beta3,PB3*length(beta3)/(max(beta3)-min(beta3)),type="l",xlab=expression(beta[3]),ylab=expression(paste(italic("p("),beta[3],"|",italic("y)"),sep="") ))

Xs<-cbind(rep(1,6),1:6,(1:6)^2)
eXB.post<- exp(t(Xs%*%t(BETAg )) )
qE<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

source("hdr2d.R")
plot.hdr2d(BETAg[,2:3],bw=c(15,15),xlab=expression(beta[2]),
           ylab=expression(beta[3]))
dev.off()


#### MH algorithm for one-sample normal problem with 

## Setup
s2<-1 
t2<-10 ; mu<-5

set.seed(1)
n<-5
y<-round(rnorm(n,10,1),2)

mu.n<-( mean(y)*n/s2 + mu/t2 )/( n/s2+1/t2) 
t2.n<-1/(n/s2+1/t2)

## MCMC
s2<-1 ; t2<-10 ; mu<-5 
y<-c(9.37, 10.18, 9.16, 11.60, 10.33)
theta<-0 ; delta<-2 ; S<-10000 ; THETA<-NULL ; set.seed(1)

for(s in 1:S)
{

  theta.star<-rnorm(1,theta,sqrt(delta))

  log.r<-( sum(dnorm(y,theta.star,sqrt(s2),log=TRUE)) +
               dnorm(theta.star,mu,sqrt(t2),log=TRUE) )  -
         ( sum(dnorm(y,theta,sqrt(s2),log=TRUE)) +
               dnorm(theta,mu,sqrt(t2),log=TRUE) ) 

  if(log(runif(1))<log.r) { theta<-theta.star }

  THETA<-c(THETA,theta)

}



#### Figure 10.3
pdf("fig10_3.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

skeep<-seq(10,S,by=10)
plot(skeep,THETA[skeep],type="l",xlab="iteration",ylab=expression(theta))

hist(THETA[-(1:50)],prob=TRUE,main="",xlab=expression(theta),ylab="density")
th<-seq(min(THETA),max(THETA),length=100)
lines(th,dnorm(th,mu.n,sqrt(t2.n)) )
dev.off()


#### MH algorithm with different proposal distributions
par(mfrow=c(4,3))
ACR<-ACF<-NULL
THETAA<-NULL
for(delta2 in 2^c(-5,-1,1,5,7) ) {
set.seed(1)
THETA<-NULL
S<-10000
theta<-0
acs<-0
delta<-2

for(s in 1:S) 
{

  theta.star<-rnorm(1,theta,sqrt(delta2))
  log.r<-sum( dnorm(y,theta.star,sqrt(s2),log=TRUE)-
              dnorm(y,theta,sqrt(s2),log=TRUE)  )  +
      dnorm(theta.star,mu,sqrt(t2),log=TRUE)-dnorm(theta,mu,sqrt(t2),log=TRUE) 

  if(log(runif(1))<log.r)  { theta<-theta.star ; acs<-acs+1 }
  THETA<-c(THETA,theta) 

}
plot(THETA[1:1000])

ACR<-c(ACR,acs/s) 
ACF<-c(ACF,acf(THETA,plot=FALSE)$acf[2]  )
THETAA<-cbind(THETAA,THETA)
}
plot(ACR,ACF) ; lines(ACR,ACF)



#### Figure 10.4 
pdf("fig10_4.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
laby<-c(expression(theta),"","","","")

for(k in c(1,3,5)) {
plot(THETAA[1:500,k],type="l",xlab="iteration",ylab=laby[k], 
    ylim=range(THETAA) )
abline(h=mu.n,lty=2)
                  }
dev.off()

THCM<-apply(THETAA,2,cumsum)
THCM<- THCM/(1:dim(THCM)[1]) 



#### Back to sparrow data
fit.mle<-glm(fledged~age+age2,family="poisson")
summary(fit.mle)

y<-fledged ; X<-cbind(rep(1,length(y)),age,age^2)
yX<-cbind(y,X)
colnames(yX)<-c("fledged","intercept","age","age2") 

n<-length(y) ; p<-dim(X)[2]

pmn.beta<-rep(0,p)
psd.beta<-rep(10,p)

var.prop<- var(log(y+1/2))*solve( t(X)%*%X )
beta<-rep(0,p)
S<-10000
BETA<-matrix(0,nrow=S,ncol=p)
ac<-0
set.seed(1)

## rmvnorm function for proposals
rmvnorm<-function(n,mu,Sigma)
{ # samples from the multivariate normal distribution
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}

## MCMC
for(s in 1:S) {

#propose a new beta

beta.p<- t(rmvnorm(1, beta, var.prop ))

lhr<- sum(dpois(y,exp(X%*%beta.p),log=T)) -
      sum(dpois(y,exp(X%*%beta),log=T)) +
      sum(dnorm(beta.p,pmn.beta,psd.beta,log=T)) -
      sum(dnorm(beta,pmn.beta,psd.beta,log=T))

if( log(runif(1))< lhr ) { beta<-beta.p ; ac<-ac+1 }

BETA[s,]<-beta
                    }
cat(ac/S,"\n")

library(coda)
apply(BETA,2,effectiveSize)



#### Figure 10.5 
pdf("fig10_5.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
blabs<-c(expression(beta[1]),expression(beta[2]),expression(beta[3]))
thin<-c(1,(1:1000)*(S/1000))
j<-3
plot(thin,BETA[thin,j],type="l",xlab="iteration",ylab=blabs[j])
abline(h=mean(BETA[,j]) )

acf(BETA[,j],ci.col="gray",xlab="lag")
acf(BETA[thin,j],xlab="lag/10",ci.col="gray")
dev.off()


#### Figure 10.6
pdf("fig10_6.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot(beta2,PB2*length(beta2)/(max(beta2)-min(beta2)) ,type="l",xlab=expression(beta[2]),ylab=expression(paste(italic("p("),beta[2],"|",italic("y)"),sep="") ) ,lwd=2,lty=2,col="gray")
lines(density(BETA[,2],adj=2),lwd=2)

plot(beta3,PB3*length(beta3)/(max(beta3)-min(beta3)),type="l",xlab=expression(beta[3]),ylab=expression(paste(italic("p("),beta[3],"|",italic("y)"),sep="") ),lwd=2,col="gray",lty=2)
lines(density(BETA[,3],adj=2),lwd=2)

Xs<-cbind(rep(1,6),1:6,(1:6)^2) 
eXB.post<- exp(t(Xs%*%t(BETA )) )
qE<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

plot( c(1,6),range(c(0,qE)),type="n",xlab="age",
   ylab="number of offspring")
lines( qE[1,],col="black",lwd=1)
lines( qE[2,],col="black",lwd=2)
lines( qE[3,],col="black",lwd=1)

dev.off()



#### Ice core example
load("icecore.RData") 

pdf("fig10_7.pdf",family="Times",height=1.75,width=5)

par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
layout(matrix( c(1,1,2),nrow=1,ncol=3) )

plot(icecore[,1],  (icecore[,3]-mean(icecore[,3]))/sd(icecore[,3]) ,
   type="l",col="black",
   xlab="year",ylab="standardized measurement",ylim=c(-2.5,3))
legend(-115000,3.2,legend=c("temp",expression(CO[2])),bty="n",
       lwd=c(2,2),col=c("black","gray"))
lines(icecore[,1],  (icecore[,2]-mean(icecore[,2]))/sd(icecore[,2]),
  type="l",col="gray")

plot(icecore[,2], icecore[,3],xlab=expression(paste(CO[2],"(ppmv)")),ylab="temperature difference (deg C)")

dev.off()



#### Figure 10.8
pdf("fig10_8.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

lmfit<-lm(icecore$tmp~icecore$co2)
hist(lmfit$res,main="",xlab="residual",ylab="frequency")
acf(lmfit$res,ci.col="gray",xlab="lag")



#### Starting values for MCMC
n<-dim(icecore)[1]
y<-icecore[,3]
X<-cbind(rep(1,n),icecore[,2])
DY<-abs(outer( (1:n),(1:n) ,"-"))

lmfit<-lm(y~-1+X)
beta<-lmfit$coef
s2<-summary(lmfit)$sigma^2
phi<-acf(lmfit$res,plot=FALSE)$acf[2]
nu0<-1 ; s20<-1 ; T0<-diag(1/1000,nrow=2)


## MCMC - 1000 scans saving every scan
set.seed(1)
S<-1000 ; odens<-S/1000
OUT<-NULL ; ac<-0 ; par(mfrow=c(1,2))
for(s in 1:S)
{

  Cor<-phi^DY  ; iCor<-solve(Cor)
  V.beta<- solve( t(X)%*%iCor%*%X/s2 + T0)
  E.beta<- V.beta%*%( t(X)%*%iCor%*%y/s2  )
  beta<-t(rmvnorm(1,E.beta,V.beta)  )

  s2<-1/rgamma(1,(nu0+n)/2,(nu0*s20+t(y-X%*%beta)%*%iCor%*%(y-X%*%beta)) /2 )

  phi.p<-abs(runif(1,phi-.1,phi+.1))
  phi.p<- min( phi.p, 2-phi.p)
  lr<- -.5*( determinant(phi.p^DY,log=TRUE)$mod -
             determinant(phi^DY,log=TRUE)$mod  +
   sum(diag( (y-X%*%beta)%*%t(y-X%*%beta)%*%(solve(phi.p^DY) -solve(phi^DY)) ) )/s2 )

  if( log(runif(1)) < lr ) { phi<-phi.p ; ac<-ac+1 }

  if(s%%odens==0)
    {
      cat(s,ac/s,beta,s2,phi,"\n") ; OUT<-rbind(OUT,c(beta,s2,phi))
    }
}

OUT.1000<-OUT
library(coda)
apply(OUT.1000,2,effectiveSize )


## MCMC - 25000 scans saving every 25th scan
set.seed(1)
S<-25000 ; odens<-S/1000
OUT<-NULL ; ac<-0 ; par(mfrow=c(1,2))
for(s in 1:S)
{

  Cor<-phi^DY  ; iCor<-solve(Cor)
  V.beta<- solve( t(X)%*%iCor%*%X/s2 + T0)
  E.beta<- V.beta%*%( t(X)%*%iCor%*%y/s2  )
  beta<-t(rmvnorm(1,E.beta,V.beta)  )

  s2<-1/rgamma(1,(nu0+n)/2,(nu0*s20+t(y-X%*%beta)%*%iCor%*%(y-X%*%beta)) /2 )

  phi.p<-abs(runif(1,phi-.1,phi+.1))
  phi.p<- min( phi.p, 2-phi.p)
  lr<- -.5*( determinant(phi.p^DY,log=TRUE)$mod -
             determinant(phi^DY,log=TRUE)$mod  +
   sum(diag( (y-X%*%beta)%*%t(y-X%*%beta)%*%(solve(phi.p^DY) -solve(phi^DY)) ) )/s2 )

  if( log(runif(1)) < lr ) { phi<-phi.p ; ac<-ac+1 }

  if(s%%odens==0)
    {
      cat(s,ac/s,beta,s2,phi,"\n") ; OUT<-rbind(OUT,c(beta,s2,phi))
    }
}

OUT.25000<-OUT
library(coda)
apply(OUT.25000,2,effectiveSize )



#### Figure 10.9
pdf("fig10_9.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(OUT.1000[,4],xlab="scan",ylab=expression(rho),type="l")
acf(OUT.1000[,4],ci.col="gray",xlab="lag")
dev.off()



#### Figure 10.10
pdf("fig10_10.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(OUT.25000[,4],xlab="scan/25",ylab=expression(rho),type="l")
acf(OUT.25000[,4],ci.col="gray",xlab="lag/25")
dev.off()



#### Figure 10.11
pdf("fig10_11.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

plot(density(OUT.25000[,2],adj=2),xlab=expression(beta[2]),
   ylab="posterior marginal density",main="")

plot(y~X[,2],xlab=expression(CO[2]),ylab="temperature")
abline(mean(OUT.25000[,1]),mean(OUT.25000[,2]),lwd=2)
abline(lmfit$coef,col="gray",lwd=2)
legend(180,2.5,legend=c("GLS estimate","OLS estimate"),bty="n",
      lwd=c(2,2),col=c("black","gray"))
dev.off()

