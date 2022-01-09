numParms = 2
startAge = 50
endAge   = 100

age = startAge:endAge
a = 50
b = 1
TT  = readRDS('Data/cCountries.rds') %>% as_tibble() %>%
  dplyr::filter(Age>=startAge,Age<=endAge,!Country%in%c('East Germany','West Germany')) %>%
  dplyr::arrange(Country,Sex,Age,Year) %>% 
  dplyr::mutate(Death = na_if(Death, 0),
                se = sqrt(1/Death),
                # Basis function
                B0 = 1,
                B1 = (Age-a)/b
  ) %>% drop_na() 

if(backtest!=FALSE) TT = TT %>% filter(Year<=backtest)

TT = TT %>% 
  group_by(Country, Sex, YoB) %>%
  filter(Sex=="Female") %>%
  dplyr::summarise(n = n()) %>% 
  filter(n==51) %>% 
  filter(YoB==min(YoB)) %>%
  dplyr::rename("minYoB"="YoB") %>%
  ungroup() %>%
  dplyr::select(-n,-Sex) %>% 
  right_join(TT) %>%
  filter(YoB>=minYoB) %>%
  dplyr::select(-minYoB)