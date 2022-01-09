source('funCohortModel.R', local = TRUE)
source('setup.R', local = TRUE)

# ==============================================================================
numSig = 1
opt = trust.optim(x0,fn=LogLikeBMAPFn,gr=LogLikeBMAPGr
                  ,TT=TT
                  ,numParms=numParms
                  ,prior=modelCountry[[1]]$prior
                  ,method="SR1"
                  ,control=list(maxit=2000,report.freq=50))
sigma = tail(opt$solution,1)
x = head(opt$solution,NROW(x0)-1)
if(randEff==1){
  #------------- Do loop settings ---------------------------
  maxAbsDis  = 10
  iter       = 0
  sigmaPrime = rep(0,2*numContr)
  sigma      = rep(sigma,2*numContr)
  while(maxAbsDis>10^(-5)){
    iter=iter+1
    cat('================================================\n')
    cat('iter         =',iter,'\n')
    cat('---------- calculate sigma ----------\n')
    opt=optim(par=sigma,fn=LogLikeSigmaFn
              ,x=x
              ,TT=TT
              ,numParms=numParms
              ,prior=modelCountry[[1]]$prior
              ,V=modelCountry[[1]]$V
              ,control=list(maxit=2000,trace=TRUE))
    sigma=opt$par
    maxAbsDis=max(abs(sigma-sigmaPrime))
    cat('converge     = ',opt$convergence, '\n')
    cat('sigma        = ',sigma, '\n')
    cat('exp(sigma)^2 = ',exp(sigma)^2, '\n')
    cat('maxAbsDis    = ',maxAbsDis, '\n')
    cat('---------- calculate x     ----------\n')
    opt=trust.optim(x,fn=LogLikeBetaFn,gr=LogLikeBetaGr
                    ,sigma=sigma
                    ,TT=TT
                    ,numParms=numParms
                    ,prior=modelCountry[[1]]$prior
                    ,V=modelCountry[[1]]$V
                    ,method="SR1"
                    ,control=list(maxit=2000,report.freq=50))
    x = opt$solution
    sigmaPrime = sigma
  }
}

z = as.matrix(c(x,sigma))
hessian  = CalHessian(x,sigma,TT,numParms,modelCountry[[1]]$prior,modelCountry[[1]]$V,randEff)
LogLike = LogLikePostFn(z,TT,numParms,modelCountry[[1]]$prior,modelCountry[[1]]$V)
AIC = 2*NROW(z) - 2*LogLike
BIC = NROW(z)*log(NROW(TT)) - 2*LogLike

result=list(x=x,
            sigma=sigma,
            hessian=hessian,
            modelCountry=modelCountry,
            LogLike=LogLike,
            AIC=AIC,
            BIC=BIC)

save(result,file=paste0("Estimates/",savePath,"/c",countries[1],sex[ss],"B0B1var",ifelse(backtest!=FALSE,backtest,""),".rData"))
