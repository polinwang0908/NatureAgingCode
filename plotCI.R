maxPeriod = 2021
p2Data = merged_all_double %>%
  dplyr::filter(Country == cc) %>%
  dplyr::select(Country, Sex, cohort, upper, lower) %>%
  dplyr::mutate(CI = "maxLife") %>%
  bind_rows(gammas %>% filter(Country == cc) %>% mutate(CI = "gamma"))

p1 = p1Data %>% filter(Country == cc) %>%
  ggplot(aes(x = B0, y = B1, color = Prior)) +
  geom_point(alpha = 0.5) +
  geom_label_repel(
    data = p1Data  %>% filter(Country == cc) %>% filter(Prior == "Prior A", cohort %in%
                                                          seq(1700, 1960, 25)),
    aes(label = cohort),
    box.padding = 0.5,
    max.overlaps = Inf,
    color = "black"
  ) +
  facet_wrap( ~ Sex, ncol = 2) +
  xlab(TeX(r'($\lambda_c$)')) +
  ylab(TeX(r'($\delta_c$)')) + theme(
    legend.position = "bottom",
    legend.title = element_text(size =
                                  9),
    legend.text = element_text(size = 9)
  )


p1



p2 = p2Data %>%
  ggplot() +
  geom_ribbon(aes(
    x = cohort,
    ymin = lower,
    ymax = upper,
    fill = CI
  ), alpha = 0.5) +
  
  facet_grid( ~ Sex) +
  geom_point(
    data = MaximumLivedPeople %>% filter(Country == cc),
    aes(x = cohort , y = maxage, color = Alive),
    alpha = 0.7
  ) +
  scale_color_manual(
    name = "Longest-lived individual per cohort",
    values = c(1, "springgreen4"),
    labels = c("Dead before 2021", "Alive in 2021"),
    drop = FALSE
  ) +
  scale_fill_manual(
    name = "95% CI’s for",
    values = c(2, 1),
    labels = list(
      TeX(r'($\Lambda_c$ (age at which mortality hazard first reaches 2/3))'),
      TeX(r'($M_c$ (maximum age at death))')
    ),
    drop = FALSE
  ) +
  scale_x_continuous(breaks = seq(1700, 1975, 20)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) +
  ylab("Years") +
  xlab("Birth year") +
  geom_polygon(
    data = data.frame(x = 1700:maxPeriod, y = (maxPeriod - 1700):0) %>% tibble::add_row(x = c(maxPeriod, 1700), y = c(Inf, Inf)),
    aes(x = x, y = y),
    fill = "grey",
    alpha = 0.4
  ) +
  coord_cartesian(xlim = c(1700, 1975), ylim = c(95, 150))

p2

p3 = p3Data %>% filter(Country == cc) %>%
  filter(cohort %in% c(seq(1700, 1970, 10), max(cohort))) %>%
  mutate(cohort = ifelse(cohort > 1960, 1970, cohort)) %>%
  pivot_longer(cols = c(compression, postponement),
               names_to = "type") %>%
  group_by(cohort, Sex, type) %>%
  dplyr::summarise(
    middle = median(value),
    upper.dist = quantile(value, 0.975, na.rm = TRUE) - middle,
    lower.dist = middle - quantile(value, 0.025, na.rm =
                                     TRUE)
  ) %>%
  group_by(Sex, cohort) %>%
  arrange(Sex, cohort, type) %>%
  dplyr::mutate(
    type = factor(type, levels = c("postponement", "compression")),
    cummean = case_when(cumprod(middle) >= 0 ~ cumsum(middle),
                        TRUE ~ middle),
    ymin = cummean - lower.dist,
    ymax = cummean + upper.dist,
    type = dplyr::recode(type, "compression" = "Due to compression", "postponement" = "Due to postponement")
  ) %>%
  ggplot() +
  geom_bar(aes(
    x = cohort,
    y = middle,
    fill = type,
    group = type
  ), stat = "identity") +
  geom_errorbar(
    aes(
      x = cohort,
      ymin = ymin,
      ymax = ymax,
      fill = type
    ),
    position = position_dodge(10),
    alpha = 0.7
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 9)) +
  
  scale_x_continuous(limits = c(1700, 1975),
                     breaks = seq(1700, 1975, 20)) +
  geom_point(
    data = CohortLE50 %>% filter(Country == cc),
    aes(x = cohort, y = `LE actual change`, shape = "Actual change in remaining life expectancy at 50")
  ) +
  facet_wrap( ~ Sex, ncol = 2) +
  xlab("Birth year (end of 10-year period)") +
  ylab(
    "Change in remaining life expectancy at age 50 \n relative to cohort born 10 years earlier (years)"
  )

p3
library(patchwork)
p = p1 + p3 + p2 + plot_annotation(title = ifelse(cc == "NewZealand", "New Zealand", cc))  + plot_layout(ncol =
                                                                                                           1)
print(p)
