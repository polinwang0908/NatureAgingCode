backtest = FALSE
source("readMortalityData.R")

TT = TT %>%
  arrange(Country, Sex, YoB, Age) %>%
  group_by(Country, Sex, YoB) %>%
  mutate(
    qx = 2 * mx / (2 + mx),
    oneMinusQx = 1 - qx,
    el = cumprod(oneMinusQx),
    el2 = case_when(Age == 50 ~ 1,
                    TRUE ~ dplyr::lag(el, n = 1, order_by = Age)),
    L = (el + el2) / 2,
    valid = all(c(50:95) %in% Age)
  ) %>% 
  dplyr::summarise(LE50 = sum(L[valid]),
                   LE50 = na_if(LE50, 0))

saveRDS(TT, file = "Data/CohortLE50.rds")
