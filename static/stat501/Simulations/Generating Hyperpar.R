
### load scenarios data
# source("multi_test_simulation_DTR.R")



### Convert lower and upper limits to alpha and beta of Gamma distribution
getalphabeta <- function(ll,ul){ 
  if(ul>1){ul <- 0.99}
  mu <- (ul+ll)/2
  # ab<- matrix(NA,nrow=50,ncol = 3)
  # for (i in 1:50) {
  #   geta <- function(a){
  #     b <- a*(1-mu)/mu
  #     abs(pbeta(ul,shape1=a,shape2=b)-pbeta(ll,shape1=a,shape2=b)-0.95)
  #   }
  #   a<- optim(par=seq(1,100,length.out=50)[i],fn=geta,method="BFGS")
  #   ab[i,] <-c(a$par,a$par*(1-mu)/mu,a$value)
  # }
  # ab<- matrix(NA,nrow=50,ncol = 3)
  # return(ab[ab[,3]==min(ab[,3]),1:2])
  
  geta <- function(a){
    b <- a*(1-mu)/mu
    abs(pbeta(ul,shape1=a,shape2=b)-pbeta(ll,shape1=a,shape2=b)-0.95)
  }
  a<- optim(par=1,fn=geta,method="BFGS")
  c(a$par,a$par*(1-mu)/mu)
  
  
}


### Convert lower and upper limits to a and b of Normal distribution
getab <- function(ll,ul){ 
  if(ul>1){ul <- 0.99}
  mu_a=(qnorm(ul)+qnorm(ll))/2
  mu_b <- 0
  sigma_a=(qnorm(ul)-qnorm(ll))/4
  LP=0.001; UP=0.999
  sigma_b<-  max(sqrt(((qnorm(UP)-mu_a)/1.96)^2 - sigma_a^2)/1.96,
                 sqrt(((mu_a-qnorm(LP))/1.96)^2 - sigma_a^2)/1.96)
  return(c(mu_a,sigma_a,mu_b,sigma_b))
}

### Explore the parameter space
ll=c(0.01,0.01,0.495,0.98) 
ul=c(0.02,0.99,0.505,0.99)
alphabeta<- matrix(NA,nrow=4,ncol =2)
ab<- matrix(NA,nrow=4,ncol =4)
for (i in 1:4) {
  alphabeta[i,] <- getalphabeta(ll[i],ul[i])
  ab[i,] <- getab(ll[i],ul[i])
}
colnames(alphabeta) <- c("alpha","beta")
rownames(alphabeta) <- c("sharp skewed right","flat center","sharp center","sharp skewed left")
colnames(ab) <- c("mu_a","sigma_a","mu_b","sigma_b")
rownames(ab) <- c("minimum mu","maximum sd","minimum sd","maximum mu")

alphabeta
ab



### Generating Hyper-parameters
genhyperpars <- function(trueinfo,
                         setting=c("pi.info","S1.info","C1.info","S2.info","C2.info"),
                         scale=c(0.8,1.2)){
  hyper.parlist <- list(fixed=c(api=2,bpi=2,
                                aS1=2,bS1=2,
                                aS2=2,bS2=2,
                                aC1=2,bC1=2,
                                aC2=2,bC2=2,
                                acovs=2,bcovs=2,
                                acovc=2,bcovc=2),
                        random=c(api=2,bpi=2,
                                 muS1=2.19,sdS1=1.16317,
                                 muS2=2.19,sdS2=1.16317,
                                 muC1=2.19,sdC1=1.16317,
                                 muC2=2.19,sdC2=1.16317,
                                 b1mu=0,b1sd=1.374,
                                 b0mu=0,b0sd=1.374)
  )
  names(trueinfo) <- c("pi","S1","S2","C1","C2")
  switch(setting,
         pi.info={
           prev <- trueinfo["pi"]
           hyper.parlist$fixed[c("api","bpi")] <- getalphabeta(ll=scale[1]*prev,ul=scale[2]*prev)
           hyper.parlist$random[c("api","bpi")] <-getalphabeta(ll=scale[1]*prev,ul=scale[2]*prev)
           return(hyper.parlist)
         },
         S1.info={
           S1 <- trueinfo["S1"]
           hyper.parlist$fixed[c("aS1","bS1")] <- getalphabeta(ll=scale[1]*S1,ul=scale[2]*S1)
           hyper.parlist$random[c("muS1","sdS1","b1mu","b1sd")]<-getab(ll=scale[1]*S1,ul=scale[2]*S1)
           return(hyper.parlist)
         },
         C1.info={
           C1 <- trueinfo["C1"]
           hyper.parlist$fixed[c("aC1","bC1")] <- getalphabeta(ll=scale[1]*C1,ul=scale[2]*C1)
           hyper.parlist$random[c("muC1","sdC1","b0mu","b0sd")]<-getab(ll=scale[1]*C1,ul=scale[2]*C1)           
           return(hyper.parlist)
         },
         S2.info={
           S2 <- trueinfo["S2"]
           hyper.parlist$fixed[c("aS2","bS2")] <- getalphabeta(ll=scale[1]*S2,ul=scale[2]*S2)
           hyper.parlist$random[c("muS2","sdS2","b1mu","b1sd")]<-getab(ll=scale[1]*S2,ul=scale[2]*S2)
           return(hyper.parlist)
         },
         C2.info={
           C2 <- trueinfo["C2"]
           hyper.parlist$fixed[c("aC2","bC2")] <- getalphabeta(ll=scale[1]*C2,ul=scale[2]*C2)
           hyper.parlist$random[c("muC2","sdC2","b0mu","b0sd")]<-getab(ll=scale[1]*C2,ul=scale[2]*C2)           
           return(hyper.parlist)
         })
}

# parnams <- c("pi","S1","S2","C1","C2")
# genhyperpars(trueinfo=sel.truth.scen[3,parnams],setting="pi.info")
# genhyperpars(trueinfo=sel.truth.scen[3,parnams],setting="S1.info")
# genhyperpars(trueinfo=sel.truth.scen[3,parnams],setting="C1.info")
# genhyperpars(trueinfo=sel.truth.scen[3,parnams],setting="S2.info")
# genhyperpars(trueinfo=sel.truth.scen[3,parnams],setting="C2.info")
