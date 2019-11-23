## ----echo=F, message=F, warning=F----------------------------------------

#get data and define prior and posterior parameters
pyg.data <- Sleuth3::ex1321
summ.pyg <- with(pyg.data,
                 by(data = Gain, INDICES = Treatment,
                    FUN = function(x){
                      data.frame(n=round(length(x),0),mean=mean(x),sd=sd(x))
                      }))

mu0 <- 0
nu <- 1
alpha <- 1/2
beta <- 100*alpha

#data
x <- pyg.data$Gain[pyg.data$Treatment=="pygmalion"]
y <- pyg.data$Gain[pyg.data$Treatment=="control"]

#data summaries
np <- summ.pyg$pygmalion$n
x.bar <- summ.pyg$pygmalion$mean

nc <- summ.pyg$control$n
y.bar <- summ.pyg$control$mean

#posterior parameters
mu.p.s <- (np/(nu+np))*x.bar
mu.c.s <- (nc/(nu+nc))*y.bar

nu.p <- nu+np
nu.c <- nu+nc

alpha.p <- alpha+np/2
alpha.c <- alpha+nc/2

beta.p <- 0.5*(sum(x^2)+2*beta-nu.p*mu.p.s^2)
beta.c <- 0.5*(sum(y^2)+2*beta-nu.c*mu.c.s^2)

S <- 10^4
rnormgamma <- function(S,a,b,mu,nu){
  #a: scale paramter for gamma component
  #b: rate paramter for gamma component
  #mu: mean parameter for normal component
  #nu: scaling factor for precision of normal component
  
  lambda <- rgamma(S,shape=a,rate=b)
  theta <- rnorm(S,mean=mu,sd=1/sqrt(nu*lambda))
  return(list(lambda=lambda,theta=theta))
}

sample.p.1 <- rnormgamma(S=S,a=alpha.p,b=beta.p,mu=mu.p.s,nu=nu.p)
sample.c.1 <- rnormgamma(S=S,a=alpha.c,b=beta.c,mu=mu.c.s,nu=nu.c)

psi <- mean((sample.p.1$theta>sample.c.1$theta)[-(1:1000)])



## ----echo=F, message=F, warning=F, fig.width=6, fig.height=4, fig.align='center'----
runnav <- cumsum(sample.p.1$theta>sample.c.1$theta)/(1:S)
ll <- runnav - 2*sqrt(psi*(1-psi)/(1:S))
ul <- runnav + 2*sqrt(psi*(1-psi)/(1:S))
plot(x=1:S,y=runnav,type="l",
     ylim=range(c(ll[-(1:100)],ul[-(1:100)])),
     xlab="number of samples (S)",
     ylab=expression(psi[S]))
abline(h=psi,col="red",lwd=2)
lines(ll,col="cornflowerblue",lty=3)
lines(ul,col="cornflowerblue",lty=3)
#legend()


## ----echo=F--------------------------------------------------------------
load("alldata")
w40 <- Y[(Y$YEAR>=1990)&(Y$FEMALE==1)&(Y$AGE==40),#&(Y$DEGREE<3),
         c("YEAR","CHILDS","AGE","DEGREE")]
w40 <- na.exclude(w40)

#empirical distribution
ct.child <- rep(0,11)
ct.child[1:length(unique(w40$CHILDS))] <- as.vector(table(w40$CHILDS))
pr.child <- ct.child/sum(ct.child)
names(pr.child) <- names(ct.child) <- 0:10

#posterior predictive NegBinom(sum(y) + a, (b+n)/(b+n+1))
a <- 2; b <- 1
n <- nrow(w40)
sum.y <- sum(w40$CHILDS)
p <- (b+n)/(b+n+1)

post.pred <- dnbinom(0:10,size=(a+sum.y),prob=p)#mu=((a+sum.y)/(b+n)))#


## ----echo=F, message=F, warning=F, fig.width=8, fig.height=4, fig.align='center'----

par(mfrow=c(1,2))
#plot empty graph device (using the option type="n")
dev <- 0.1
plot( x=c(-dev,10+dev), y=c(0,.4),xlab="number of children",
      ylab=expression(p(Y[i]==y[i])),type="n",main="")

#add barplots for theta=0.05, 0.1 and 0.2
points(0:10-dev,pr.child,type="h",col="black",lwd=4)
points(0:10+dev,post.pred,type="h",col=gray(.75),lwd=4)
legend(4,0.35,legend=c("empirical",
                       "predictive"),
       bty="n", lwd=c(2,2),col=c("black",gray(.75)))

S <- 10000
O21 <- rep(NA,S)
theta.k <- ytilde.ks <- NULL

for(k in 1:S){
  theta.k <- rgamma(1,a+sum.y,b+n)
  ytilde.ks <- rpois(n,theta.k)
  O21[k] <- sum(ytilde.ks==2)/sum(ytilde.ks==1)
}
hist(x=O21,main="",xlab="posterior predictive odds for 2 vs 1 children",freq=F)
abline(v=pr.child[3]/pr.child[2],lwd=3,col="red")


## ----echo=F, message=F, warning=F, fig.width=6, fig.height=4, fig.align='center'----
S <- 5000
psi <- pcauchy(2,lower.tail = F)
th.sd <- sqrt(psi*(1-psi)/(1:S))

#direct MC
x.k <- rcauchy(S)
ra.x <- cumsum(x.k>2)/(1:S)
ll.x <- ra.x - 2*th.sd
ul.x <- ra.x + 2*th.sd

#using uniforms MC
y.k <- runif(S,0,2)
z.k <- 2/(pi*(1+y.k^2))
sd.z <- sqrt(0.0285/(1:S))
ra.y <- 0.5-cumsum(z.k)/(1:S)
ll.y <- ra.y - 2*sd.z
ul.y <- ra.y + 2*sd.z

plot(x=1:S,y=ra.x,type="l",
     ylim=c(0.1,0.2),
     xlab="number of samples (S)",
     ylab=expression(psi[S]),
     col="blue",lwd=2)
lines(ll.x,col="blue",lty=2)
lines(ul.x,col="blue",lty=2)

lines(ra.y,col="orange",lwd=2)
lines(ll.y,col="orange",lty=2)
lines(ul.y,col="orange",lty=2)
abline(h=psi,col="black",lwd=2)
legend("topright",bty="n",legend=c("true","MC Cauchy","MC Unif"),col=c("black","blue","orange"),lty=1)


## ----echo=F, message=F, warning=F, fig.width=8, fig.height=4, fig.align='center'----
S <- 10000
Y <- runif(S,0,10)
gY <- 10*exp(-2*abs(Y-5))
g <- function(x){10*exp(-2*abs(x-5))}

par(mfrow=c(1,2))
plot(density(gY,from=0, to=10),type="l",
      main="",xlab=expression(z),ylab=expression(p(z)))
curve(expr=g,from=0.001,to=10,xlab="y",ylab="g(y)")


## ----echo=F, message=F, warning=F, fig.width=5, fig.height=4, fig.align='center'----
g <- function(x) 10*exp(-2*abs(x-5))

dispvals <- c(0.1,0.5,1,2)
yvals <- seq(0,10,by=0.1)
plot(yvals,dcauchy(yvals,location = 5,scale=0.1)/3,lty=2,col=2,
     type="n",ylab="rescaled g(y) and q(y)",xlab="y")
lines(yvals,g(yvals)/10,type="l",lwd=2)
for(dd in dispvals){
  lines(yvals,dcauchy(yvals,location = 5,scale=dd)/3,lty=1,col="cornflowerblue",lwd=1)
  lines(yvals,dnorm(yvals,mean = 5,sd=dd)/3,lty=1,col="orange",lwd=1)
}
legend("topright",lty=rep(1,2),col=c("cornflowerblue","orange","black"),legend=c(expression(paste("Cauchy(",5,",",gamma,")")),expression(paste("Normal(",5,",",sigma^2,")")),"g(y)"),bty="n")


## ----echo=TRUE,results="asis",warning=FALSE------------------------------
IS.fn <- function(X,dq){
  g(X)*dunif(X,0,10)/dq(X)
} 

S <- 10000

results.normal <- results.cauchy <-list()
results.unif <- g(runif(S,0,10))
uu.res <- c(estimate=mean(results.unif),var=var(results.unif))

k <- 1
for(dd in dispvals){
  dqn <- function(X){dnorm(X,mean=5,sd=dd)}
  dqc <- function(X){dcauchy(X,location=5,scale=dd)}
  
  results.normal[[k]] <- IS.fn(X=rnorm(S,mean=5,sd=dd),dq=dqn)
  results.cauchy[[k]] <- IS.fn(X=rcauchy(S,location=5,scale=dd),dq=dqc)
  k <- k+1
}

results.all <- data.frame(rbind(uu.res,
                     do.call(rbind,lapply(results.normal,
                                          function(xx)c(estimate=mean(xx),var=var(xx)))),
                     do.call(rbind,lapply(results.cauchy,
                                          function(xx)c(estimate=mean(xx),var=var(xx))))))
results.all <- cbind(family=c("Uniform",rep("Normal",4),rep("Cauchy",4)),
                     disp=c("",rep(dispvals,2)),
                     results.all)

  
library(xtable)
print(xtable(results.all,
             caption="Importance sampling estimate and estimator variance 
             for different proposal distributions.",
             label="tab:ISestvar"),
      include.rownames = FALSE,booktabs = T)


## ----echo=F, message=F, warning=F, fig.width=8, fig.height=4, fig.align='center'----
S <- 10000
Y <- runif(S,0,5)
hY <- 10*exp(-2*abs(Y-5))
g <- function(x){10*exp(-2*abs(x-5))}

par(mfrow=c(1,2))
plot(density(hY,from=0, to=10),type="l",
      main="",xlab=expression(z),ylab=expression(p(z)))
curve(expr=g,from=0.001,to=5,xlab="y",ylab="h(y)")


## ----echo=TRUE,results="asis",warning=FALSE------------------------------
g.U <- function(x){10*exp(-2*abs(x-5))}

S <- 10000

IS.weight <- function(X,q){
  dunif(X,0,5)/q(X)
}


X.unif <- runif(S,0,5)
g.unif <- g.U(X.unif)

mu <- 5; sigma <- 1
qnormal <- function(x){dnorm(x,mean=mu,sd=sigma)}
X.normal <- rnorm(S,mean=mu,sd=sigma)
g.normal <- g.U(X.normal)*IS.weight(X.normal,qnormal)

Y.half <-  X.normal
Y.half[Y.half>5] <- 10-Y.half[Y.half>5]
qhalf <- function(x){2*dnorm(x,mean=mu,sd=sigma)}
g.half <- g.U(Y.half)*IS.weight(Y.half,qhalf)

mu <- 5; sigma <- 1
qcauchy <- function(x){dcauchy(x,location=mu,scale=sigma)}
X.cauchy <- rcauchy(S,location = mu, scale=sigma)
g.cauchy <- g.U(X.cauchy)*IS.weight(X.cauchy,qcauchy)

Y.half2 <-  X.cauchy
Y.half2[Y.half2>5] <- 10-Y.half2[Y.half2>5]
qhalf2 <- function(x){2*dcauchy(x,location=mu)}
g.half2 <- g.U(Y.half2)*IS.weight(Y.half2,qhalf2)


results.all2 <- round(rbind(MC.unif = c(estimate=mean(g.unif),var=var(g.unif)),
      IS.nor = c(estimate=mean(g.normal),var=var(g.normal)),
      IS.half = c(estimate=mean(g.half),var=var(g.half)),
      IS.cauchy = c(estimate=mean(g.cauchy),var=var(g.cauchy)),
      IS.half2 = c(estimate=mean(g.half2),var=var(g.half2))),4)
results.all2 <- cbind(family=c("Uniform","Normal","Half.Normal",
                               "Cauchy","Half.Cauchy"),
                      disp=c(NA,1,1,1,1),location=c(NA,rep(5,4)),
                      results.all2)

library(xtable)
print(xtable(results.all2,
             caption="Importance sampling estimate and estimator variance
             for different proposal distributions.",
             label="tab:ISestvar2"),
      include.rownames = FALSE,booktabs = T)

