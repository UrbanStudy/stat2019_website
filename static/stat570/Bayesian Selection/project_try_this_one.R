#this code generally works. you can change the predictors you want to test 
#in line 25 to try different combinations of predictors in the model.


source("Bayesian Selection/VarSelectHC.R")
source("Bayesian Selection/summaryout.R")

all.data <- read.csv("Willamette_habitat_features.csv",header=T)
all.data[is.na(all.data)]<-0
ind.insample <- sample(1:178,100)
X <- data.frame(all.data[-c(1,2)])
#not sure if these operations make sense, but based on our conversation with James 
#it seems that area_2m nests area_1m, and similarly for perc_2m and perc_1m (this last I'm not so sure)
X$area_2m <- X$area_2m-X$area_1m
X$perc_2m <- X$perc_2m-X$perc_1m

#for the response, we can log or sqrt the response so that it lives in the whole real line
y <- all.data[2]
orig.namvars <- names(X)
names(X) <- paste0("x",1:ncol(X))

#model blows up with all predictors probably due multicolinearity, so I just 
#sample 10 predictors to test. Determine which predictors are strongly/perfectly colinear and 
#get rid of some of them
# vtest <- c("x3","x5","x6","x7","x9","x12","x15")
vtest <- c("x5","x6","x7","x12")
datain <- data.frame(y=y[ind.insample,],X[ind.insample,vtest])
data.holdout <- data.frame(y=y[-ind.insample,],X[-ind.insample,vtest])
modpriorvec=c("EPP","HOP","HTP")


theformula <- as.formula(paste("y~",paste0(vtest,collapse="+")))

#notice that if I specify max.deg=2, in full.fomula I only need to specify the main effects,
#and max.deg=2 takes care of inotrducing all main effects and interactions.

res=VarSelectHC(full.formula=theformula,
                data=datain,
                base.formula=as.formula(. ~ 1),
                maxdeg=2,
                nodes.to.remove=NULL,
                model.prior.type=modpriorvec,
                model.prior.pars = "children",
                beta.prior.type = "IP",
                beta.prior.pars = list(alpha=1,nu=1),
                niter=5000)

summary.res <- summaryout(mcmc.out=res,insampledata=datain,modelprior.nams=modpriorvec,
                          shr.adj=T,outsampledata=data.holdout,respnam="y",top.ave=10,betaprtype="IP",
                          parsprbeta=list(alpha=1,nu=1))
