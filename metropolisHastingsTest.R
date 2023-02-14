metropolisHastingsTest = function(dd,ff='var',numSim=10^5,report=10^3,backtest=FALSE){
  source("setCountries.R")
  Countries = countries
  Sex = c('Male','Female')
  ccc = dd%%NROW(Countries)
  cc = countries = Countries[ifelse(ccc==0,NROW(Countries),ccc)]
  ss = Sex[(dd-1)%/%NROW(Countries)+1]
  print(paste0("Running ",cc, ss))
  #####################################
  # MH
  #####################################
  sortT<-function(TT) data.table(TT[order(TT$Age,TT$Year),])
  CalFitted<-function(beta0,TT,numParms){
    beta    = data.table(t(matrix(beta0,numParms)))
    cohort  = unique(TT$YoB[order(TT$YoB)])
    beta    = cbind(cohort,beta)
    colnames(beta)=c('YoB',paste('beta',(0:(numParms-1)),sep=''))
    TT     = merge(TT,beta,by='YoB')
    TT      = sortT(TT)
    beta    = t(Matrix(as.matrix(subset(TT,select=paste('beta',(0:(numParms-1)),sep='')))))
    B       = t(Matrix(as.matrix(subset(TT,select=paste('B'   ,(0:(numParms-1)),sep='')))))
    fitted  = colSums(B*beta)
    return(fitted)
  }
  parmsTransToOptim = function(x0, numSig, fn = "mean"){
    # fn == "mean" or "sum"
    temp = split(x0, rep(countries, numParms*numAllCohort+numSig))
    xMatrix = list()
    sigma = list()
    for(dd in 1:numContr){
      sigma[[dd]] = tail(temp[[dd]],numSig)
      xMatrix[[dd]] = matrix(head(temp[[dd]],NROW(temp[[dd]])-numSig),numParms)
    }
    if(fn == "mean")slope = mean(do.call(cbind,xMatrix)[2,])
    if(fn == "sum") slope = sum(do.call(cbind,xMatrix)[2,])
    out = c(do.call(c,lapply(1:numContr, function(dd)c(xMatrix[[dd]][-2,],sigma[[dd]]))),slope)
  }
  
  parmsTransBack = function(x0, numSig, type = "list", sep = FALSE){
    # sep -- saparate x and sigma
    temp = list()
    x = list()
    sigma = list()
    for(dd in 1:numContr){
      temp[[dd]] = x0[1:((numParms-1)*modelCountry[[dd]]$number+numSig)]
      x0 = x0[((numParms-1)*modelCountry[[dd]]$number+numSig+1):NROW(x0)]
      xMatrix = matrix(head(temp[[dd]],NROW(temp[[dd]])-numSig),numParms-1)
      xMatrix = rbind(xMatrix[1,],
                      rep(x0[NROW(x0)],modelCountry[[dd]]$number),
                      xMatrix[(numParms-1)*(numParms!=2),])
      x[[dd]] = c(xMatrix)
      sigma[[dd]] = tail(temp[[dd]],numSig)
      temp[[dd]] = c(x[[dd]], sigma[[dd]])
    }
    if(type!="list")temp = do.call(c,temp)
    if(type!="list"&sep)temp = list(x=do.call(c,x),sigma=do.call(c,sigma))
    if(type=="list"&sep)temp = lapply(1:numContr,function(dd)list(x=x[[dd]],sigma=sigma[[dd]]))
    temp
  }
  LogLikeFn<-function(data){
    z = data
    x = head(z, NROW(z)-2)
    sigma = tail(z, 2)
    cohort = unique(TT$YoB[order(TT$YoB)])
    fitted = CalFitted(x,TT,numParms)
    f      = log(TT$mx)-fitted
    
    sigmavec = c(exp(sigma[1])^2+TT$se^2)
    missing  = is.na(sigmavec)
    covM     = Diagonal(x=sigmavec)+V*exp(sigma[2])^2
    f =f [!missing]
    covM=covM[!missing,!missing]
    
    L=chol(covM)
    logdetsigma=2*sum(log(diag(L)))
    delta.x=prior$gamma%*%x-prior$mu
    yy=0.5*logdetsigma+0.5*t(f)%*%solve(covM,f)+0.5*t(delta.x)%*%solve(prior$sigma,delta.x)
    return(yy=as.numeric(yy))
  }
  #################################
  library(tidyverse)
  library(Matrix)
  library(data.table)
  dump   = 1
  source("readMortalityData.R")
  
  TT = TT %>% filter(Country==cc,Sex==ss)
  
  load(paste0('Estimates/',startAge,endAge,'/priorA/c',cc,ss,'B0B1',ff,ifelse(backtest!=FALSE,backtest,''),'.rData'))
  x = result$x
  sigma = result$sigma
  hessian = result$hessian
  modelCountry = result$modelCountry
  prior = modelCountry[[1]]$prior
  V = Matrix(modelCountry[[1]]$V)
  numAllCohort = modelCountry[[1]]$number
  numContr = 1
  
  if(ff=='fix') x = parmsTransToOptim(x, numSig = 0, fn = "mean")
  
  invhess = solve(hessian)
  posdefinvshess = (invhess + t(invhess)) / 2
  sigmat = posdefinvshess
  if(ff=='fix') {
    z = as.matrix(c(x[-NROW(x)],sigma,x[NROW(x)]))
  } else {
    z = as.matrix(c(x,sigma))
    }
  accepted = 0
  jump = 2.4/NROW(z)^.5
  zSample = matrix(NA,NROW(z),numSim/dump)
  j = 0
  
  for(i in (j+1):(numSim)){
    z_c = MASS::mvrnorm(1,z[,NCOL(z)],jump^2*sigmat)
    if(ff=='fix'){
      r = exp(-LogLikeFn(parmsTransBack(z_c,numSig=2,type="not"))
              +LogLikeFn(parmsTransBack(z[,NCOL(z)],numSig=2,type="not")))    
    } else {
      r = exp(-LogLikeFn(z_c)
              +LogLikeFn(z[,NCOL(z)]))    
    }
    
    if (runif(1)<r){
      z = cbind(z,z_c)
      accepted = accepted + 1
    } else {
      z = cbind(z,z[,NCOL(z)])
    }
    z = z[,(NCOL(z)-1):NCOL(z)]
    if(i%%dump==0)zSample[,i] = z[,NCOL(z)]
    if(i%%report==0){
      show(paste0(cc,ss,i,'\n'))
      show(paste0('acceptance = ',accepted/i,'\n'))
      result = list(zSample=zSample,accepted=accepted)
      save(result,file=paste0('MH/',startAge,endAge,'MHc',cc,ss,ff,ifelse(backtest!=FALSE,backtest,''),'.rData'))
    }
  }
  
  return(NULL)
}
