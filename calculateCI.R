library(haven)
library(tidyverse)
library(optimr)
library(data.table)
source("sampleMH.R")

if(backtest!=FALSE){
  countries = "Sweden"
} else {
  if(robust==TRUE){
    countries = "USA"
  }
  source("setCountries.R")
}

temp = read_dta("Data/Population50.dta") %>%
  mutate(cohort = Year - Age) %>%
  pivot_longer(
    cols = c(Female1, Male1),
    names_to = "Sex",
    values_to = "Population"
  ) %>%
  mutate(Sex = dplyr::recode(Sex, "Female1" = "Female", "Male1" = "Male")) %>%
  dplyr::select(cohort, Country, Sex, Population) %>%
  right_join(map_dfr(c(countries),function(cc){
    map_dfr(c("Male","Female"),function(ss){
      load(paste0('Estimates/',startAge,endAge,'/priorA/c',cc,ss,'B0B1var',ifelse(backtest!=FALSE,backtest,''),'.rData'))
      numParms = 2
      modelCountry = result$modelCountry
      load(paste0('MH/',startAge,endAge,'/MHc',cc,ss,'var',ifelse(backtest!=FALSE,backtest,''),'.rData'))
      z = result$zSample
      Total = sum(!is.na(z[1,]))
      sample = z[,seq(burnin + 1,Total,everyTh)]
      xSample = map_dfr(1:dim(sample)[2],function(i){
        parms = t(matrix(sample[1:(NROW(sample)-2),i],2)) %>%
          data.frame() %>% 
          `colnames<-`(c("B0","B1")) %>%
          mutate(Country = cc,
                 Sex = ss,
                 cohort = modelCountry[[1]]$cohort,
                 sim = i
          )
      })    
    })
  }) ) %>%
  data.table()

temp = temp %>% 
  mutate(p = exp((exp(B0)-mumax)/B1),
         stdp = sqrt(Population*p*(1-p)),
         maxint = floor(Population*p+5*stdp),
         Xi = 50+(log(mumax)-B0)/B1)

options(dplyr.summarise.inform = FALSE)

maxLifeCDF = function(x, data){
  p = data %>% 
    rowwise() %>%
    dplyr::mutate(temp = case_when(x<Xi~0,
                                   TRUE~log(1-exp(-mumax*(x-Xi)))),
                  points = list(dbinom(0:maxint,floor(Population),p)),
                  interim = list(case_when(temp==0~0,
                                           TRUE~exp(seq(0,temp*maxint,temp)))),
                  maxLifeCDF = (sum(points*interim)-points[1]*interim[1])/(sum(points)-points[1])) %>%
    group_by() %>%
    dplyr::summarise(maxLifeCDF = sum(maxLifeCDF)/n()) %>%
    pull(maxLifeCDF)
  p
}

optimCDF = function(x, q, data){
  abs(q - maxLifeCDF(x, data))*10^10
  
}

getPercentile = function(q, data){
  if(is.na(data$Population[1])){
    out = NA
  } else {
    this.x = mean(data$Xi)
    upper = 150
    lower = 80
    segments = 3
    
    error = 10
    iter = 0
    while(error>0.001&iter<=7){
      x = seq(lower,upper,(upper-lower)/segments)
      fValues = sapply(x,FUN=maxLifeCDF, data=data)
      
      upper = min(x[which(cumsum(fValues>q)==1)])
      lower = min(x[which(cumsum(fValues>q)==1)-1])
      this.x = (upper + lower)/2
      error = abs(maxLifeCDF(this.x, data) - q)
      iter = iter + 1
    }
    
    out = this.x
  }
}

for(cc in countries){
  for(ss in c("Male","Female")){
    print(paste0("calculating ",cc, " ",ss))
    try(temp %>% 
        dplyr::filter(Country==cc,Sex==ss) %>%
        dplyr::group_by(Country,Sex,cohort) %>%
          tidyr::nest()  %>% 
        dplyr::mutate(lower = map_dbl(data, ~ getPercentile(0.025, data=.x)),
                      upper = map_dbl(data, ~ getPercentile(0.975, data=.x))) %>%
        dplyr::select(-data) %>%
          saveRDS(file=paste0(paste0("Estimates/",startAge,endAge,"/CI/",cc,ss,"CI",ifelse(backtest!=FALSE,backtest,''),ifelse(mumax==1,'One',''),".rds")))    )
  }
}
