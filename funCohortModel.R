LogLikeFn <- function(beta0,TT,numParms){
  beta     = Matrix(beta0[-length(beta0)],numParms,NROW(TT))
  sigma    = exp(beta0[length(beta0)])
  x        = Matrix(as.matrix(subset(TT,select=paste('B',(0:(numParms-1)),sep=''))))
  fitted   = colSums(t(x)*beta)
  f        = log(TT$mx)-fitted
  intsigma = sqrt(sigma^2+TT$se^2)
  
  y        = sum(-log(intsigma)-.5*(f/intsigma)^2,na.rm = T)
  y        = -y
  return(y)
}
LogLikeGr <- function(beta0,TT,numParms){
  beta     = Matrix(beta0[-length(beta0)],numParms,NROW(TT))
  sigma    = exp(beta0[length(beta0)])
  x        = Matrix(as.matrix(subset(TT,select=paste('B',(0:(numParms-1)),sep=''))))
  fitted   = colSums(t(x)*beta)
  f        = log(TT$mx)-fitted
  intsigma = sqrt(sigma^2+TT$se^2)
  
  xx            = t(Matrix(f/intsigma^2,nrow(TT),numParms))
  xx[is.na(xx)] = 0
  d             = rep(0,numParms+1)
  d[1:numParms] = diag(xx%*%x)
  d[numParms+1] = sum(-(sigma/intsigma)^2+(sigma*f/intsigma^2)^2,na.rm=T)
  return(d=-d)
}
LogLikeBMAPFn <- function(beta0,TT,numParms,prior){
  cohort   = unique(TT$YoB[order(TT$YoB)])
  sigma    = beta0[length(beta0)]
  beta     = t(Matrix(beta0[-length(beta0)],numParms))
  yy       = sum(sapply(1:nrow(beta),function(r)LogLikeFn(c(beta[r,],sigma),TT[TT$YoB==cohort[r],],numParms)))
  delta.x  = prior$gamma%*%beta0[-length(beta0)]-prior$mu
  
  yy       = yy+as.numeric(0.5*t(delta.x)%*%solve(prior$sigma,delta.x))
  return(yy)
}
LogLikeBMAPAllFn = function(beta0,TT,numParms,modelCountry){
  # fix slope -> matrix(intercept, curvartuare, slope, sigma)
  beta0 = parmsTransBack(beta0,numSig)
  do.call(sum,llply(1:numContr,function(dd){
    LogLikeBMAPFn(beta0 = beta0[[dd]],TT = TT %>% filter(Country == countries[dd]),numParms = numParms, prior = modelCountry[[dd]]$prior)
  }, .parallel = FALSE))
}
LogLikeBMAPGr <- function(beta0,TT,numParms,prior){
  cohort   = unique(TT$YoB[order(TT$YoB)])
  sigma    = beta0[length(beta0)]
  beta     = t(Matrix(beta0[-length(beta0)],numParms))
  gg      = t(sapply(1:nrow(beta),function(r)LogLikeGr(c(beta[r,],sigma),TT[TT$YoB==cohort[r],],numParms)))
  delta.x = prior$gamma%*%beta0[-length(beta0)]-prior$mu
  
  gg=cbind(matrix(t(gg[,1:numParms]),1)+t(t(prior$gamma)%*%solve(prior$sigma,delta.x)),sum(gg[,numParms+1]))
  return(as.numeric(gg))
}
LogLikeBMAPAllGr = function(beta0,TT,numParms,modelCountry){
  beta0 = parmsTransBack(beta0,numSig)
  temp = do.call(c,llply(1:numContr,function(dd){
    LogLikeBMAPGr(beta0 = beta0[[dd]],TT = TT %>% filter(Country == countries[dd]),numParms = numParms, prior = modelCountry[[dd]]$prior)
  }, .parallel = FALSE))
  parmsTransToOptim(temp, numSig, fn = "sum")
}
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
LogLikeBetaFn<-function(x,sigma,TT,numParms,prior,V){
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
LogLikeBetaAllFn = function(x,sigma,TT,numParms,modelCountry){
  x = parmsTransBack(x,numSig=0,type="list",sep=TRUE)
  sigma   = split(sigma, rep(countries, 2))
  do.call(sum,lapply(1:numContr,function(dd){
    LogLikeBetaFn(x = x[[dd]]$x, 
                  sigma = sigma[[dd]],
                  TT = TT %>% filter(Country == countries[dd]),
                  numParms = numParms, 
                  prior = modelCountry[[dd]]$prior,
                  V = modelCountry[[dd]]$V)
  }))
}
TransformParms<-function(x)  cbind(t(matrix(x[-length(x)],numParms)),rep(x[length(x)],(length(x)-1)/numParms))
LogLikeBetaGr<-function(x,sigma,TT,numParms,prior,V){
  cohort = unique(TT$YoB[order(TT$YoB)])
  fitted = CalFitted(x,TT,numParms)
  f      = log(TT$mx)-fitted
  cohort = unique(TT$YoB[order(TT$YoB)])
  B      = subset(TT,select=paste('B',(0:(numParms-1)),sep=''))
  
  sigmavec = c(exp(sigma[1])^2+TT$se^2)
  missing  = is.na(sigmavec)
  covM     = Diagonal(x=sigmavec)+V*exp(sigma[2])^2
  f    = f[!missing]
  TT   = TT[!missing,]
  B    = B[!missing,]
  covM = covM[!missing,!missing]
  
  delta.x=prior$gamma%*%x-prior$mu
  F = solve(covM,f)
  gg = c(sapply(1:length(cohort),function(i)colSums(B[TT$YoB==cohort[i],,drop=FALSE]*F[TT$YoB==cohort[i],drop=FALSE],na.rm=TRUE)))
  dd=-gg+t(prior$gamma)%*%solve(prior$sigma,delta.x)
  return(as.numeric(dd))
}
LogLikeBetaAllGr = function(x,sigma,TT,numParms,modelCountry){
  x = parmsTransBack(x,numSig=0,type="list",sep=TRUE)
  sigma   = split(sigma, rep(countries, 2))
  temp = do.call(c,llply(1:numContr,function(dd){
    LogLikeBetaGr(x = x[[dd]]$x,
                  sigma = sigma[[dd]],
                  TT = TT %>% filter(Country == countries[dd]),
                  numParms = numParms, 
                  prior = modelCountry[[dd]]$prior,
                  V = modelCountry[[dd]]$V)
  }, .parallel = FALSE))
  parmsTransToOptim(temp, numSig = 0, fn = "sum")
}
LogLikeSigmaFn<-function(sigma,x,TT,numParms,prior,V){
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
LogLikeSigmaAllFn = function(sigma,x,TT,numParms,modelCountry){
  x = parmsTransBack(x,numSig=0,type="list",sep=TRUE)
  sigma   = split(sigma, rep(countries, 2))
  do.call(sum,llply(1:numContr,function(dd){
    LogLikeSigmaFn(sigma = sigma[[dd]],
                   x = x[[dd]]$x, 
                  TT = TT %>% filter(Country == countries[dd]),
                  numParms = numParms, 
                  prior = modelCountry[[dd]]$prior,
                  V = modelCountry[[dd]]$V)
  }, .parallel = FALSE))
}
CalHessian<-function(x,sigma,TT,numParms,prior,V,randEff){
  library(sparseHessianFD)
  bdlength = ifelse(randEff,length(x),length(x)+1)
  bd  = Matrix(TRUE,bdlength,bdlength)
  mc  = Matrix.to.Coord(tril(bd))
  if(randEff==1){
    obj = sparseHessianFD(x,fn=LogLikeBetaFn,gr=LogLikeBetaGr
                          ,rows=mc$rows
                          ,cols=mc$cols
                          ,sigma=sigma
                          ,TT=TT
                          ,numParms=numParms
                          ,prior=prior
                          ,V=V
                          #,fixSlope=fixSlope
                          )  
    hess = obj$hessian(x)
    epsilon = 1e-5
    sigma1Plus  = sigma[1] + .5*epsilon
    sigma2Plus  = sigma[2] + .5*epsilon
    sigma1Minus = sigma[1] - .5*epsilon
    sigma2Minus = sigma[2] - .5*epsilon
    AA=rbind( LogLikeBetaGr(sigma=c(sigma1Plus ,sigma[2]   ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)
             -LogLikeBetaGr(sigma=c(sigma1Minus,sigma[2]   ),x=x,TT=TT,numParms=numParms,prior=prior,V=V),
              LogLikeBetaGr(sigma=c(sigma[1]   ,sigma2Plus ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)
             -LogLikeBetaGr(sigma=c(sigma[1]   ,sigma2Minus),x=x,TT=TT,numParms=numParms,prior=prior,V=V))/epsilon
    
    B11   =   LogLikeBetaFn(sigma=c(sigma1Plus ,sigma[2]   ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)-
            2*LogLikeBetaFn(sigma=c(sigma[1]   ,sigma[2]   ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)+
              LogLikeBetaFn(sigma=c(sigma1Minus,sigma[2]   ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)
    
    B12=B21 = LogLikeBetaFn(sigma=c(sigma1Plus ,sigma2Plus ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)-
              LogLikeBetaFn(sigma=c(sigma1Plus ,sigma2Minus),x=x,TT=TT,numParms=numParms,prior=prior,V=V)-
              LogLikeBetaFn(sigma=c(sigma1Minus,sigma2Plus ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)+
              LogLikeBetaFn(sigma=c(sigma1Minus,sigma2Minus),x=x,TT=TT,numParms=numParms,prior=prior,V=V)
    
    B22     = LogLikeBetaFn(sigma=c(sigma[1]   ,sigma2Plus ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)-
            2*LogLikeBetaFn(sigma=c(sigma[1]   ,sigma[2]   ),x=x,TT=TT,numParms=numParms,prior=prior,V=V)+
              LogLikeBetaFn(sigma=c(sigma[1]   ,sigma2Minus),x=x,TT=TT,numParms=numParms,prior=prior,V=V)
    
    BB=rbind(c(B11,B12),c(B21,B22))/(.5*epsilon)^2
    hess=rbind(cbind(hess,t(AA)),cbind(AA,BB))
  }
  if(randEff==0){
    obj = sparseHessianFD(c(x,sigma),fn=LogLikeBMAPFn,gr=LogLikeBMAPGr
                          ,rows=mc$rows
                          ,cols=mc$cols
                          ,TT=TT
                          ,numParms=numParms
                          ,prior=prior
                          
                          )  
    hess = obj$hessian(c(x,sigma))
    
  }
  #hess=nearPD(rbind(cbind(hess,t(AA)),cbind(AA,BB)))
  #hess=hess$mat
  return(hess)
}
LogLikePostFn<-function(z,TT,numParms,prior,V){
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
  yy=0.5*logdetsigma+0.5*t(f)%*%solve(covM,f)#+0.5*t(delta.x)%*%solve(prior$sigma,delta.x)
  return(yy=-as.numeric(yy))
}
sortT<-function(TT) data.table(TT[order(TT$Age,TT$Year),])
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
