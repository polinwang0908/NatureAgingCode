randEff = 1
sex = c('Male','Female')
# =======================================================================================================
# ---------------------------------------- Read data ----------------------------------------------------
source("readMortalityData.R")
savePath = case_when(pp=="A"~"priorA",
                     pp=="B"~"priorB",
                     TRUE~"priorC")

TT = TT %>% filter(Sex==sex[ss])
countries = unique(TT$Country)
numContr = NROW(countries)

modelCountry = lapply(1:numContr, function(i) {
  TT = TT %>% filter(Country == countries[i])
  #------- Set up age random effects -------------
  V = TT %>% lm(log(mx)~as.factor(Age), data = .) %>% sparse.model.matrix()
  V = cbind(rowSums(V)==1,V[,-1])
  V = V%*%t(V)
  
  period = min(TT$Year):max(TT$Year)
  cohort = min(TT$YoB):max(TT$YoB)
  number = NROW(cohort)
  if(pp=="A")cohortPrior = min(TT$YoB):(max(TT$Year)-90)  
  if(pp!="A") cohortPrior = max(1869,min(TT$YoB)):(max(TT$Year)-90)
  coByCo = ts(t(sapply(cohort,function(i){
    tryCatch(lm(log(mx)~B1,data=TT[TT$YoB==i,])$coef, error = function(e) rep(NA,numParms))
  })
  ),min(cohort))
  
  # --------------------------------------------------------------------------------
  # VAR
  startC = min(cohortPrior)
  endC   = max(cohortPrior)
  cc     = which(cohort==startC):which(cohort==endC)
  
  parms      = window(coByCo,startC,endC)
  varData    = diff(parms)
  
  result = list(period = period,
                cohort = cohort,
                cohortPrior = cohortPrior,
                coByCo = coByCo,
                number = number,
                varData = varData,
                V = V)
})

if(pp=='C'){
  varDataAll = map_dfr(1:NROW(modelCountry),function(i){
    modelCountry[[i]]$varData %>% data.frame() %>% `colnames<-`(c("B0","B1"))
  })}

modelCountry = lapply(1:numContr, function(i) {
  cohort = modelCountry[[i]]$cohort
  number = modelCountry[[i]]$number
  cohortPrior = modelCountry[[i]]$cohortPrior
  coByCo = modelCountry[[i]]$coByCo
  
  if(pp=='C') varData = varDataAll
  if(pp!='C') varData = modelCountry[[i]]$varData
  
  x0.mini    = rep(0,NROW(varData))
  varModel   = vars::VAR(varData)#,exogen=x0.mini)
  kappa.mini = rep(0,numParms)
  gamma.mini = vars::Bcoef(varModel)[,1:numParms]
  mu.mini    = vars::Bcoef(varModel)[,numParms+1]
  sigma.mini = ar.ols(varData,order=1)$var.pred
  
  L = cbind(Diagonal((number-1)*numParms,-1),matrix(0,(number-1)*numParms,numParms))+
    cbind(matrix(0,(number-1)*numParms,numParms),Diagonal((number-1)*numParms,1))
  
  prior = list()
  prior$mu    = rep(mu.mini,number-2)
  prior$sigma = bdiag(replicate(number-2, sigma.mini,simplify=FALSE))
  prior$gamma = (cbind(bdiag(replicate(number-2,-gamma.mini,simplify=FALSE)),matrix(0,(number-2)*numParms,numParms))+
                   cbind(matrix(0,(number-2)*numParms,numParms),Diagonal((number-2)*numParms)))%*%L
  
  startC = min(cohortPrior)
  endC   = max(cohortPrior)
  cc     = which(cohort==startC):which(cohort==endC)
  
  x0 = t(replicate(NROW(cohort),apply(coByCo,2,mean,na.rm=TRUE),simplify = TRUE))
  x0[cc,] = window(coByCo,startC,endC)
  x0 = c(t(x0),1)
  
  result = modelCountry[[i]]
  result$prior = prior
  result$varData = varData
  result$x0 = x0
  
  result
})

modelCountry = lapply(1,function(i)modelCountry[[ifelse(backtest,which(countries=="Sweden"),dd)]])
TT = TT %>% filter(Country==ifelse(backtest,"Sweden",countries[dd]))
numContr = 1
countries = countries[ifelse(backtest,which(countries=="Sweden"),dd)]

x0 = do.call(c,lapply(1:numContr,function(dd)modelCountry[[1]]$x0))
numAllCohort = sapply(1:numContr,function(dd)modelCountry[[1]]$number)