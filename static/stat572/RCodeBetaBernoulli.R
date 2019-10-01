#--------------------------------------------------------------
# Prior and posterior from the Beta-Bernoulli model (coin tossing)
#--------------------------------------------------------------

#--------------------------------------------------------------
#plotting and playing with the Beta(a,b) prior distribution
#--------------------------------------------------------------
genprBetas <- function(a,b,colset="blue",addF=F){
  curve(dbeta(x,shape1=a,shape2=b),from=0.0,to=1,
        ylab=expression(p(theta)),xlab=expression(theta),
        n=500,col=colset,add=addF,lwd=2,ylim=c(0,5))
  abline(v=0.51,lty=3,lwd=1)
}


genprBetas(a=1,b=1,colset="grey",addF=F)
genprBetas(a=2,b=2,colset="blue",addF=T)
genprBetas(a=10,b=10,colset="red",addF=T)
genprBetas(a=0.5,b=0.5,colset="darkgreen",addF=T)
genprBetas(a=0.5,b=2,colset="magenta",addF=T)
genprBetas(a=2,b=0.5,colset="cornflowerblue",addF=T)
# abline(v=0.51,lwd=1,lty=3)
legend("top",lty=rep(1,4),lwd=2,
       col=c("grey","blue","red","darkgreen","magenta","cornflowerblue"),
       legend=c("a=1,b=1","a=2,b=2","a=10,b=10","a=0.5,b=0.5","a=0.5,b=2","a=2,b=0.5"),
       ncol=2,cex=0.7)


#--------------------------------------------------------------
#plotting and playing with the Beta(a+sum(y),b+n-sum(y)) posterior distribution
#--------------------------------------------------------------
genBetas <- function(n,colset="blue",addF=F){
  y<-rbinom(n,1,0.51)
  curve(dbeta(x,shape1=(1+sum(y)),shape2=(1+n-sum(y))),from=0.0,to=1,
        ylab=expression(p(paste(theta,"|",y[1:n]))),n=500,col=colset,add=addF,lwd=2)
  abline(v=0.51,lty=3,lwd=1)  
}


genBetas(1000,colset="cornflowerblue",addF=F)
genBetas(100,colset="red",addF=T)
genBetas(10,colset="darkgreen",addF=T)
genBetas(0,colset="black",addF=T)
abline(v=0.51,lwd=1,lty=3)
legend("topright",lty=rep(1,4),lwd=2,
       col=c("cornflowerblue","red","darkgreen","black"),
       legend=c("n=1,000","n=100","n=10","n=0"))