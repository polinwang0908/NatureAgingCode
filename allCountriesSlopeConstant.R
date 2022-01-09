source('funCohortModelSlopeConstant.R', local = TRUE)
source('setup.R', local = TRUE)
  
x0 = parmsTransToOptim(x0,1,fn="mean")
# =========================================================================

numSig = 1
opt = trust.optim(x0,fn=LogLikeBMAPAllFn,gr=LogLikeBMAPAllGr
                ,TT=TT
                ,numParms=numParms
                ,modelCountry=modelCountry
                ,method="SR1"
                ,control=list(maxit=2000,report.freq=50))
x     = parmsTransBack(opt$solution,numSig=1,type="not",sep=TRUE)
sigma = x$sigma
x     = parmsTransToOptim(x$x, numSig = 0, fn = "mean")
if(randEff==1){
  #------------- Do loop settings ---------------------------
  numSig = 2
  maxAbsDis  = 10
  iter       = 0
  sigmaPrime = rep(0,2*numContr)
  sigma      = rep(sigma,2*numContr)
  while(maxAbsDis>10^(-5)){
    iter=iter+1
    cat('================================================\n')
    cat('iter         =',iter,'\n')
    cat('---------- calculate sigma ----------\n')
    opt=optim(par=sigma,fn=LogLikeSigmaAllFn
               ,x=x
               ,TT=TT
               ,numParms=numParms
               ,modelCountry=modelCountry
               ,control=list(maxit=2000,trace=TRUE))
    sigma=opt$par
    maxAbsDis=max(abs(sigma-sigmaPrime))
    cat('converge     = ',opt$convergence, '\n')
    cat('sigma        = ',sigma, '\n')
    cat('exp(sigma)^2 = ',exp(sigma)^2, '\n')
    cat('maxAbsDis    = ',maxAbsDis, '\n')
    cat('---------- calculate x     ----------\n')
    opt=trust.optim(x,fn=LogLikeBetaAllFn,gr=LogLikeBetaAllGr
                    ,sigma=sigma
                    ,TT=TT
                    ,numParms=numParms
                    ,modelCountry=modelCountry
                    ,method="SR1"
                    ,control=list(maxit=2000,report.freq=50))
    x = opt$solution
    sigmaPrime = sigma
  }
}
hessian  = CalHessian(x,sigma,TT,numParms,modelCountry,randEff)
z = as.matrix(c(x[-NROW(x)],sigma,x[NROW(x)]))
LogLike = LogLikePostFn(parmsTransBack(z,numSig=2,type="not"),TT,numParms,modelCountry[[1]]$prior,modelCountry[[1]]$V)
AIC = 2*NROW(z) - 2*LogLike
BIC = NROW(z)*log(NROW(TT)) - 2*LogLike


x = parmsTransBack(x,numSig=0,type="list",sep=TRUE)
x = x[[1]]$x

result = list(
  x = x,
  sigma = sigma,
  modelCountry = modelCountry,
  hessian = hessian,
  LogLike = LogLike,
  AIC = AIC,
  BIC = BIC
)
save(result,file=paste0("Estimates/",savePath,"/c",countries[1],sex[ss],"B0B1fix.rData"))
