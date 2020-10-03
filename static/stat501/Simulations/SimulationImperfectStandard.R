library(ggplot2)
library(ggthemes)
library(tidyverse)
library(purrr)
library(viridis)
library(rjags)
library(parallel)
library(foreach)
library(doParallel)
library(furrr)

load("DataSims.RData")
source("Generating Hyperpar.R")

keep.scenarios <- sim[,"scenario"]%in%sel.truth.scen[,"id"]
sims_sel <- sim[keep.scenarios,]
cvars <- c("pt11","pt10","pt01","pt00")
parnams <- c("pi","S1","S2","C1","C2")

makedatalist <- function(type=c("fixed","random"),listcomps){
  switch(type,
         fixed={
           data.fixed <- vector(mode = "list",length=(length(listcomps[[3]])+2))
           nampars <- names(listcomps[[3]])
           names(data.fixed)=c("t12","N",nampars)
           data.fixed[[1]] <- listcomps[[1]]
           data.fixed[[2]] <- listcomps[[2]]
           for(k in nampars) data.fixed[[k]] <- unname(listcomps[[3]][k])         
           data.fixed
         },
         random={
           data.rnd <- vector(mode = "list",length=(length(listcomps[[4]])+3))
           nampars <- names(listcomps[[4]])
           names(data.rnd)=c("n4","res1", "res2",nampars)
           data.rnd[[1]] <- listcomps[[1]]
           data.rnd[[2]] <- listcomps[[2]]
           data.rnd[[3]] <- listcomps[[3]]
           for(k in nampars) data.rnd[[k]] <- unname(listcomps[[4]][k])       
           data.rnd
         })
}

run_allmethods <- function(idrun){
  
  #scenario number
  scen.num <- unname(sims_sel[idrun,"scenario"])
  
  #true parameter info for the current scenario
  sceninfo <- sel.truth.scen[sel.truth.scen[,"id"]==scen.num, parnams]
  
  #data generated for particular run of scenario considered
  rundata <- unname(sims_sel[idrun,c("T11","T10","T01","T00")]) 
  
  #get hyperpar lists for all settings
  hyper.settings <- c("pi.info","S1.info","C1.info")#,"S2.info","C2.info")
  hyperpars.all <- hyper.settings %>% map(genhyperpars,trueinfo=sceninfo)
  names(hyperpars.all) <- hyper.settings
  
  
  future::plan(multiprocess)
  hyperpars.all %>% future_map(function(x){
    N <- sum(rundata)
    data.fixed <- makedatalist(type="fixed",
                               listcomps=list(t12=rundata,N=N,
                                              hyperpars=x$fixed))
    
    #fixed model
    m.fixed <- jags.model("jagsCode/fixedJAGS.jags", data = data.fixed)
    m.fixed.samps <- coda.samples(m.fixed, c("s1", "s2", "c1", "c2", "prev",
                                             "covs12","covc12"),n.iter=20000)  
    #random model
    data.random <- makedatalist(type="random",
                                listcomps=list(n4=N,
                                               res1 = rep(c(1,1,0,0),times=rundata),
                                               res2 = rep(c(1,0,1,0),times=rundata),
                                               hyperpars=x$random))
    
    m.random <- jags.model("jagsCode/randomJAGS.jags", data = data.random)
    m.random.samps <- coda.samples(m.random, c("s1", "s2","c1","c2","prev"), n.iter=20000)
    
    return(list(res.fixed=summary(m.fixed.samps),
                res.random=summary(m.random.samps)))
    
  }, .progress = TRUE)
}



#forking parallelization, only works on mac or linux machines!!!!

# numCores <- detectCores()
results_allsims <- lapply(1:nrow(sims_sel), run_allmethods)
# mclapply(1:nrow(sims_sel), run_allmethods, mc.cores = numCores)

# #socket parallelization, works on any type of machine
# registerDoParallel(numCores)  # use multicore, set to the number of our cores
# foreach (i=1:nrow(sims_sel)) %dopar% { run_allmethods(i) }
