library(ggplot2)
library(ggthemes)
library(tidyverse)
library(purrr)
library(RANN)
# Each scenario contains 'nsim' times of simulations
# Each simulation include a vector of "Y" with 0-1, another vector of "T" with 1-4
# The summary "sim" contains 48*nsim rows and 6 colomns of "Y0","Y1","T11","T10","T01","T00"


N <- 400
nsim <- 5

cuts <- 8
nscen <- prod(rep(cuts,5))
Ylist  <- Tlist <- lapply(1:nscen,function(a)matrix(NA,ncol=nsim,nrow=N))

pivec <- seq(0.05,0.9,length.out = cuts)
S1vec <- seq(0.5,0.9,length.out = cuts)#c(0.6,0.8)
C1vec <- seq(0.5,0.9,length.out = cuts)#c(0.6,0.8)
S2vec <- seq(0.5,0.9,length.out = cuts)#c(0.6,0.8)
C2vec <- seq(0.5,0.9,length.out = cuts)#
kappa.s <- kappa.c <- NA
scenario <- matrix(NA,ncol=13,nrow=length(Ylist))

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
            Tlist[[i]][,h] <- sample(1:4,N,prob = scenario[i,6:9], replace=T)
          }
          i=i+1          
        }  
      }  
    } 
  }  
}
colnames(scenario) <- c("pi","S1","S2","C1","C2","pt11","pt10","pt01","pt00","py11","py10","py01","py00")
# scenario
#str(Ylist) 
#str(labellist)
calctable <- function(vv){
  c("0"=sum(vv==0),"1"=sum(vv==1),
    "2"=sum(vv==2),"3"=sum(vv==3),
    "4"=sum(vv==4))
}
sim <- do.call(rbind,lapply(1:length(Ylist),function(ss){
  mm <- cbind(rep(ss,nsim),apply(Ylist[[ss]],2,table) %>% t(.),
              t(apply(Ylist[[ss]]*Tlist[[ss]],2,calctable)[-1,]),
              t(apply(Tlist[[ss]],2,calctable)[-1,]))
  colnames(mm) <- c("scenario","Y0","Y1","Y11","Y10","Y01","Y00","T11","T10","T01","T00")
  mm
  }))

  


Aphasia<- readxl::read_xlsx("binary_vectors.xlsx")
testname <- c("paralg","human")
subjectname <- c("Lexicality","Phonology","Semantic")
datalist <- list(
  Lexicality = table(Aphasia$Lexicality_paralg ,Aphasia$Lexicality_human,dnn=testname),
  Phonology = table(Aphasia$Phonology_paralg ,Aphasia$Phonology_human,dnn=testname),
  Semantic = table(Aphasia$Semantic_paralg ,Aphasia$Semantic_human,dnn=testname)
)

realdata <- do.call(rbind,lapply(datalist,
                                 function(ss){res <- as.vector(ss/sum(ss))
                                 names(res)=c("T00","T01","T10","T11")
                                 res}))
#T1: human, T2:paralg
realdata <- realdata[,c("T11","T10","T01","T00")]
cvars <- c("pt11","pt10","pt01","pt00")
scenario <- cbind(scenario,id=1:nrow(scenario)) 
nn <- 1
sel.truth.scen <- scenario[nn2(scenario[,cvars],realdata,k=nn)$nn.idx,]

namtype = rownames(realdata)
toplot.scen <- data.frame(type=rep(namtype,nn),
                          rep=rep(1:nn,each=3),id=1:nrow(sel.truth.scen),
                          sel.truth.scen[,cvars],row.names=NULL) %>%
  gather(key=var,value=probs,-type,-id,-rep)

toplot.true <- data.frame(type=namtype,realdata,row.names = NULL) %>% 
  gather(key=var,value=probs,-type) %>% 
  mutate(var=paste0("p",tolower(var)))

ggplot(toplot.scen,aes(var,probs)) + 
  geom_line(aes(group=id)) +
  geom_line(data=toplot.true,colour="firebrick",aes(group=type)) +
  facet_wrap(~type)

sim <- sim[sim[,"scenario"]%in%sel.truth.scen[,"id"],]
#prune number of scenarios by removing scenarios that yield similar info:
# cvars <- c("pt11","pt10","pt01","pt00")
# set.seed <- 123445
# fit.cluster0 <- kmeans(scenario[,cvars], 9) # 4 cluster solution
# scenario <- cbind(scenario,cluster=fit.cluster0$cluster,id=1:nrow(scenario)) 
# 
# clust.centers <- cbind(cluster=1:9,fit.cluster0$centers) %>% as_tibble %>% 
#   gather(var,prob,-cluster)
# 
# scenario[,c(cvars,"cluster","id")] %>% 
#   as_tibble() %>% 
#   gather(var,prob,-cluster,-id) %>% 
#   ggplot(aes(var,prob)) + 
#   geom_line(aes(color=id,group=id)) +
#   geom_line(data=clust.centers,colour="firebrick",aes(group=cluster)) +
#   facet_wrap(~cluster)
# 
# 
# 
# 
# clust.centers <- fit.cluster0$centers
# sel.truth.scen <- scenario[nn2(scenario[,cvars],clust.centers,k=1)$nn.idx,]
# #keep one scenario from those that yield similar info (by symmetry): 
# #i.e., 4&5, 7&8. Split clusters that have sufficiently diff behavior: 9
# sel.truth.scen <- sel.truth.scen[-c(8,9),]
# # sel.truth.scen <- rbind(sel.truth.scen,
# #                         scenario[scenario[,"cluster"]==9,][1:2,])
# sel.truth.scen[,"cluster"] <- 1:nrow(sel.truth.scen)

save(list=c("sel.truth.scen","sim"),file="DataSims.RData")
