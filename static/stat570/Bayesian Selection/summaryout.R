summaryout <-
function(mcmc.out,insampledata,modelprior.nams,shr.adj=T,
                       outsampledata,respnam,top.ave=10,betaprtype="IP",
                       parsprbeta=list(alpha=1,nu=1)){
  shrinkagemaker<-function(type,n,p0,eff=TRUE,w.pars=list(alpha=1,nu=1)){
    if(type=="IP"){
      function(Rsq,p){
        q<-p*eff+1
        cf<-function(x){
          x^3*(1-x)*(
            (n-p)*q*(1-x^3)/(n+q*(1-x^3)^2)-
              (n-p0)*q*(1-x^3)/((1-Rsq)*n+q*(1-x^3)^2)-
              1/(2*(2-x^3))
          )+
            (
              (p-p0)*x^3/(1+x+x^2)+
                (5*x-3)/(6)
            )
        }
        sol<-uniroot(cf,c(0,1),tol=1e-20)
        a<-1/sol$root-1
        am<-1-sol$root^3
        am2<-am^(2)
        m<-(n-p)/2*log(n+q*am2)-(n-p0)/2*log((1-Rsq)*n+q*am2)+(p-p0)/2*log(am2)+1/6*log(1-am)-1/2*log(1+am)+2*log(1+(a-1)*(1-am)^(1/3))
        integrand<-function(t){
          u<-(1-t)/(1-t+t*a)
          v<-1-u^3
          w<-v^2
          exp((n-p)/2*log(n+q*w)-(n-p0)/2*log((1-Rsq)*n+q*w)+(p-p0)/2*log(w)+1/2*log(u)-1/2*log(1+v)+2*log(1+(a-1)*u)-m)
        }
        integral<-integrate(integrand,0,1,subdivisions=1000L,rel.tol=1E-10)
        logbf<-log(integral$value)+log(6/pi)-log(a)+m+(p-p0)/2*log(q)
        integrand<-function(t){
          u<-(1-t)/(1-t+t*a)
          v<-1-u^3
          w<-v^2
          exp(log(n)-log(n+w*q)+(n-p)/2*log(n+q*w)-(n-p0)/2*log((1-Rsq)*n+q*w)+(p-p0)/2*log(w)+1/2*log(u)-1/2*log(1+v)+2*log(1+(a-1)*u)-m)
        }
        integral<-integrate(integrand,0,1,subdivisions=1000L,rel.tol=1E-10)
        logint<-log(integral$value)+log(6/pi)-log(a)+m+(p-p0)/2*log(q)
        exp(logint-logbf)
      }
    }else if(type=="NP"){
      if(!(w.pars$alpha>0)){return("alpha<=0 for the NP")}
      alpha<-w.pars$alpha
      function(Rsq,p){
        #-(p-p0)/2*log(n/((1+p*eff)*alpha)+1)-(n-p0)/2*log(1-Rsq*(1/(1+alpha*(1+p*eff)/n)))
        exp(log(n)-log(n+alpha*(1+p*eff)))
      }  
    }else if(type=="ZS"){
      if(!(w.pars$alpha>0)){return("alpha<=0 for the ZS")}
      alpha<-w.pars$alpha
      function(Rsq,p){
        q<-p*eff+1
        cpar<-c(
          -alpha*q^2,
          q^2-alpha*q*n*(2-Rsq),
          q*n*(2-Rsq)-alpha*n^2*(1-Rsq)-n*q*(Rsq*(n-p0)-(p-p0)),
          n^2*(1-Rsq)+n^2*(p-p0)*(1-Rsq)
        )
        cpar2<-c(
          cpar[4]-cpar[3]+cpar[2]-cpar[1],
          cpar[3]-2*cpar[2]+3*cpar[1],
          cpar[2]-3*cpar[1],
          cpar[1]
        )
        cf<-function(x){
          cpar2[1]*x^3+cpar2[2]*x^2+cpar2[3]*x+cpar2[4]
        }
        sol<-uniroot(cf,c(0,1),tol=1e-20)
        a<-(1/sol$root-1)^(-.5)
        am2<-a^(-2)
        m<-(n-p)/2*log(n+q*am2)-(n-p0)/2*log((1-Rsq)*n+q*am2)+(p-p0)/2*log(am2)-alpha*am2/2+2*log(2)
        integrand<-function(t){
          v<-1/a*(1/t-1)
          w<-v^2
          exp((n-p)/2*log(n+q*w)-(n-p0)/2*log((1-Rsq)*n+q*w)+(p-p0)/2*log(w)-alpha*w/2+2*log(1+a/(1/v+(t==0)))-m)
        }
        integral<-integrate(integrand,0,1,subdivisions=1000L,rel.tol=1E-10)
        logbf<-log(integral$value)+0.5*log(2*alpha/pi)-log(a)+m+(p-p0)/2*log(q)
        integrand<-function(t){
          v<-1/a*(1/t-1)
          w<-v^2
          exp(log(n)-log(n+w*q)+(n-p)/2*log(n+q*w)-(n-p0)/2*log((1-Rsq)*n+q*w)+(p-p0)/2*log(w)-alpha*w/2+2*log(1+a/(1/v+(t==0)))-m)
        }
        integral<-integrate(integrand,0,1,subdivisions=1000L,rel.tol=1E-10)
        logint<-log(integral$value)+0.5*log(2*alpha/pi)-log(a)+m+(p-p0)/2*log(q)
        exp(logint-logbf)
      }
    }else if(type=="HG"){
      if(!(w.pars$alpha>0)){return("alpha<=0 for the HG")}
      if(!(w.pars$nu>0)){return("nu<=0 for the HG")}
      alpha<-w.pars$alpha
      nu<-w.pars$nu
      function(Rsq,p){
        alpha=w.pars$alpha
        nu=w.pars$nu
        q<-p*eff+1
        
        b<-1/nu+1
        tb<-2^b
        tbm1<-tb-1
        bp1ob<-(b+1)/b
        
        cpar<-c(
          tbm1^7*bp1ob*alpha*q^2+tb*tbm1^6*(-alpha*(nu+1)*q^2),
          tbm1^5*bp1ob*(q^2+alpha*q*n*(2-Rsq))+
            tb*tbm1^4*(-(nu+1)*alpha*q*n*(2-Rsq)-alpha*n*q*(Rsq*(n-p0)-(p-p0))),
          tbm1^3*bp1ob*(q*n*(2-Rsq)+alpha*n^2*(1-Rsq))+
            tb*tbm1^2*(-n*q*(Rsq*(n-p0)-(p-p0))-alpha*(1+nu)*n^2*(1-Rsq)+alpha*n^2*(p-p0)*(1-Rsq)),
          tbm1*bp1ob*n^2*(1-Rsq)+tb*n^2*(p-p0)*(1-Rsq)
        )
        cpar2<-c(
          cpar[4]-cpar[3]+cpar[2]-cpar[1],
          cpar[3]-2*cpar[2]+3*cpar[1],
          cpar[2]-3*cpar[1],
          cpar[1]
        )
        cf<-function(x){
          cpar2[1]*x^3+cpar2[2]*x^2+cpar2[3]*x+cpar2[4]
        }
        sol<-uniroot(cf,c(0,1),tol=1e-20)
        a<-(1/sol$root-1)^(-.5)
        am2<-a^(-2)*tbm1^2
        m<-(n-p)/2*log(n+q*am2)-(n-p0)/2*log((1-Rsq)*n+q*am2)+(p-p0)/2*log(am2)-(nu+1)/2*log(1+alpha*am2)+(1+b)*log(2)
        integrand<-function(t){
          v<-1/a*((1-t)^(-b)-1)
          w<-v^2
          exp((n-p)/2*log(n+q*w)-(n-p0)/2*log((1-Rsq)*n+q*w)+(p-p0)/2*log(w)-(nu+1)/2*log(1+alpha*w)+(1+b)/b*log(1+a/(1/v+sign(t-1)+1))-m)
        }
        integral<-integrate(integrand,0,1,subdivisions=1000L,rel.tol=1E-10)
        logbf<-log(integral$value)+0.5*log(alpha)+log(2)+log(b)-log(a)-lbeta(0.5,0.5*nu)+m+(p-p0)/2*log(q)
        integrand<-function(t){
          v<-1/a*((1-t)^(-b)-1)
          w<-v^2
          exp(log(n)-log(n+w*q)+(n-p)/2*log(n+q*w)-(n-p0)/2*log((1-Rsq)*n+q*w)+(p-p0)/2*log(w)-(nu+1)/2*log(1+alpha*w)+(1+b)/b*log(1+a/(1/v+sign(t-1)+1))-m)
        }
        integral<-integrate(integrand,0,1,subdivisions=1000L,rel.tol=1E-10)
        logint<-log(integral$value)+0.5*log(alpha)+log(2)+log(b)-log(a)-lbeta(0.5,0.5*nu)+m+(p-p0)/2*log(q)
        exp(logint-logbf)
      }
    }  
  }
  
  nodes <- mcmc.out$nodes
  p0 = length(nodes[[1]]$nodes)
  shparfn <- shrinkagemaker(type=betaprtype,n=nrow(insampledata),
                            p0=p0,eff=TRUE,w.pars=parsprbeta)
  colResp <- which(colnames(insampledata)==respnam)
  colNum <- which((unlist(lapply(insampledata,class))!="factor"))
  colNum <- colNum[colNum!=colResp]
  
  insampledata[,colNum] = scale(insampledata[,colNum])
  outsampledata[,colNum] = scale(outsampledata[,colNum])

  
  #get posteriors
  postall <- do.call(rbind,lapply(mcmc.out$models,function(mod){c(mod$modpr+mod$logBF,recursive=T)}))
  rownames(postall)=NULL
  postall <- apply(postall,2,function(a)exp(a-max(a)))
  postall <- apply(postall,2,function(b)b/sum(b))
  
  #--- extract top models for each prior, make model formula, fit, predict and model average RMSPE
  posTopmods <- apply(postall,2,function(x)order(x,decreasing=T)[1:top.ave])
  postTopmods <- do.call(cbind,lapply(1:ncol(postall),function(x)postall[posTopmods[,x],x]))
  cumpostTop <- apply(postTopmods,2,sum)
  Topmods <- apply(posTopmods,2,function(x)mcmc.out$models[x])
  colnames(postTopmods) = modelprior.nams
  
  formTopmods  <- lapply(Topmods,function(pr){
    lapply(pr,function(mod){
      if(length(mod$nodes)!=1){
        preds = paste0(unlist(lapply(nodes[mod$nodes[-1]],"[","formula")),
                       collapse="+")  
      }else{
        preds = "1"
      }
      
      as.formula(paste(paste(respnam,"~"),preds))  
    })
  })
  
  namesTopmods  <- lapply(Topmods,function(pr){
    nam=unlist(lapply(pr,function(mod){
      if(length(mod$nodes)!=1){
        paste0(unlist(lapply(nodes[mod$nodes[-1]],"[","name")),
               collapse="+")  
      }else{
        "1"
      }
      }))
    names(nam)=NULL
    return(nam)
    })
  
  listTopModels <- lapply(1:ncol(postTopmods),function(pr){
    forms=namesTopmods[[pr]]
    data.frame(Mod=forms,Post=postTopmods[,pr])
  })
  names(listTopModels) = modelprior.nams
  MSPE.ave <- unlist(lapply(1:length(formTopmods),function(prnum){
    pr <- formTopmods[[prnum]]
    postprobs <- postTopmods[,prnum]/sum(postTopmods[,prnum])
    w.resids<- rowSums(do.call(cbind,lapply(1:top.ave,function(modnum){
      form <- pr[[modnum]]
      fitmod<-lm(formula=form,data=data.frame(insampledata))
      pp <- length(fitmod$coefficients)
      if(shr.adj& pp>1){
        shr.fact <- shparfn(Rsq=summary(fitmod)$r.squared,p=pp)
        Xr <- model.matrix(terms(form),insampledata)[,-1]#as.matrix(fitmod$model[,-1])
        nn <- (fitmod$df+pp)
        b0 = (sum((1/nn)*Xr%*%matrix(( (1-shr.fact)*fitmod$coefficients[-1]),ncol=1))+
                fitmod$coefficients[1])
        fitmod$coefficients = c(b0,shr.fact*fitmod$coefficients[-1])
      }
      return((outsampledata[,respnam]-
                predict(fitmod,newdata=data.frame(outsampledata),type="response",se.fit=T)$fit)*postprobs[modnum])
      
    })))
    return(sqrt(mean(w.resids^2)))
  }))
  names(MSPE.ave) <- modelprior.nams
  
  
  #--- Extract HPM for each prior, make formula, fit and predict
  postHPM <- apply(postall,2,max)
  HPM <- lapply(1:ncol(postall),function(x)mcmc.out$models[[which(postall[,x]==max(postall[,x]))]])
  namesHPMs <- unlist(lapply(HPM,function(mod){
    if(length(mod$nodes)>1){
      paste0(c(lapply(nodes[mod$nodes[-1]],"[","name"),recursive=T),
             collapse=",")
    }else{
      "1"
    }
  }))
  
  names(namesHPMs) = names(HPM) = names(cumpostTop) = modelprior.nams
  
  formHPMs  <- lapply(HPM,function(mod){
    if(length(mod$nodes)>1){
      as.formula(paste(paste(respnam,"~"),paste0(unlist(lapply(nodes[mod$nodes[-1]],"[","formula")),
                                                 collapse="+")))
    }else{
      as.formula(paste(respnam,"~ 1"))
    }
  })
  
  fittedHPM <- lapply(formHPMs,FUN=function(ff){
    fitmod = lm(ff,data=data.frame(insampledata))
    pp <- length(fitmod$coefficients)
    if(shr.adj& pp>1){
      shr.fact <- shparfn(Rsq=summary(fitmod)$r.squared,p=pp)
      Xr <- model.matrix(terms(ff),insampledata)[,-1]#as.matrix(fitmod$model[,-1])
      nn <- (fitmod$df+pp)
      b0 = (sum((1/nn)*Xr%*%matrix(( (1-shr.fact)*fitmod$coefficients[-1]),ncol=1))+
              fitmod$coefficients[1])
      fitmod$coefficients = c(b0,shr.fact*fitmod$coefficients[-1])
    }
    return(fitmod)
  })
  
  pred.HPM <- sapply(fittedHPM,FUN=predict,newdata=data.frame(outsampledata),type="response",se.fit=T)
  PredValues=data.frame(pred.HPM[1,],check.names=F)
  colnames(PredValues)=modelprior.nams
  
  MSPE.HPM=apply(PredValues,2,FUN=function(x) sqrt(mean((outsampledata[,respnam]-x)^2)))
  Rsq <- unlist(lapply(HPM,"[","Rsq"))
  names(Rsq)=modelprior.nams
  
  return(list(formulaHPMs=namesHPMs,TopModels=listTopModels,
              post.HPM=postHPM,postcumm.Top=cumpostTop,MSPE.HPM=MSPE.HPM,
              MSPE.ave=MSPE.ave,Rsq.HPMs=Rsq))
}
