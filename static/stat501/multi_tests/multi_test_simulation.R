# Generate 48 scenario with all the combination of levels of pi, S1,S2,C1,C2
# Each scenario contains 'nsim' times of simulations
# Each simulation include a vector of "Y" with 0-1, another vector of "T" with 1-4
# The summary "sim" contains 48*nsim rows and 6 colomns of "Y0","Y1","T11","T10","T01","T00"


N <- 1000
nsim <- 1
Ylist  <- Tlist <- lapply(1:48,function(a)matrix(NA,ncol=nsim,nrow=N))

pivec <- c(0.2,0.5,0.8)
S1vec <- c(0.2,0.8)
C1vec <- c(0.2,0.8)
S2vec <- c(0.2,0.8)
C2vec <- c(0.2,0.8)
kappa.s <- kappa.c <- NA
scenario <- matrix(NA,ncol=13,nrow=48)

calculateprobs <- function(pi,S1,S2,C1,C2){
  kappa.s <- (min(S1,S2)-S1*S2)/2
  kappa.c <- (min(C1,C2)-C1*C2)/2  
  
  py11 <- (S1*S2+kappa.s)
  py10 <- (S1*(1-S2)-kappa.s)
  py01 <- ((1-S1)*S2-kappa.s)
  py00 <- ((1-S1)*(1-S2)+kappa.s)
  
  pt11 <- pi*(S1*S2+kappa.s)+(1-pi)*((1-C1)*(1-C2)+kappa.c)
  pt10 <- pi*(S1*(1-S2)-kappa.s)+(1-pi)*((1-C1)*C2-kappa.c)
  pt01 <- pi*((1-S1)*S2-kappa.s)+(1-pi)*(C1*(1-C2)-kappa.c)
  pt00 <- pi*((1-S1)*(1-S2)+kappa.s)+(1-pi)*(C1*C2+kappa.c)
  return(c(pt11,pt10,pt01,pt00,py11,py10,py01,py00))
}

n.pivec <- length(pivec)
Y4list1 <- lapply(1:(48/n.pivec),function(a)matrix(NA,ncol=nsim,nrow=N*pivec[1]))
Y4list2 <- lapply(1:(48/n.pivec),function(a)matrix(NA,ncol=nsim,nrow=N*pivec[2]))
Y4list3 <- lapply(1:(48/n.pivec),function(a)matrix(NA,ncol=nsim,nrow=N*pivec[3]))
Y4list <- c(Y4list1,Y4list2,Y4list3)
# str(Y4list)

i <-1
for(j in seq_along(pivec)){
  for(k in seq_along(S1vec)){
    for(l in seq_along(C1vec)){
      for(m in seq_along(S2vec)){
        for(n in seq_along(C2vec)){
          scenario[i,1:5] <- c(pivec[j],S1vec[k],C1vec[l],S2vec[m],C2vec[n])
          scenario[i,6:13] <-  calculateprobs(pivec[j],S1vec[k],C1vec[l],S2vec[m],C2vec[n])    
          for(h in 1:nsim){
            Ylist[[i]][,h] <- rbinom(N,1,pivec[j])
            Tlabel <- rmultinom(N, size = 1,prob = scenario[i,6:9])
            Ylabel <- rmultinom(N*pivec[j], size = 1,prob = scenario[i,10:13])
            for (t in 1:4){            
              Tlist[[i]][,h][which(Tlabel[t,]==1)]<-t # c("11","10","01","00")[t]
              Y4list[[i]][,h][which(Ylabel[t,]==1)]<-t
            }
          }
          i=i+1          
        }  
      }  
    } 
  }  
}
colnames(scenario) <- c("pi","S1","S2","C1","C2","pt11","pt10","pt01","pt00","py11","py10","py01","py00")
scenario
#str(Ylist) 
#str(labellist)


sim <- matrix(NA,ncol=10,nrow=(48*nsim))
s <- 1
for(i in 1:48){
  for(h in 1:nsim){ 
    sim[s,1:2] <- table(Ylist[[i]][,h])
    sim[s,3:6] <- table(Y4list[[i]][,h])[1:4]
    sim[s,7:10] <- table(Tlist[[i]][,h])[1:4]
    s <- s+1
  }
}
colnames(sim) <- c("Y0","Y1","Y11","Y10","Y01","Y00","T11","T10","T01","T00")
sim
