
# Script to calculate the probability that IXCHIQs benefits are greater than its risks. 
# File outputs Figures 4 and S6, S7, S8.  

library(Hmisc)
library(truncnorm)
library(tidyverse)
library(cowplot)
library(patchwork)

dir.create("output", showWarnings = FALSE)
source("R/functions.R")


cols = c("#FFA500", "#FF1493", "#756BB1","#17A2B8")

my_theme = theme_bw() +
  theme(
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.spacing.y = unit(0, "cm"),  # Reduced to zero
    legend.spacing.x = unit(0, "cm"),  # Added to control horizontal spacing
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "right")

all_risk = read.csv("data/all_risk.csv")[,-1]
all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) 
n = 100000  # posterior samples
pop = 10000 # cases averted per pop vaccinated 

# define baseline parameters

VE = 0.95

# endemic_lambda = 0.007 # Kang paper 
endemic_lambda = 0.024 # Gabriel's paper
endemic_AR = (1-(1-endemic_lambda)^10)
low_AR = 0.05
high_AR = 0.3
days_hol = 10
# traveller in Largee outbreak
# 0.3 = 1- exp(-09 * r)
r = - log(0.7)/90
traveller_AR = 1  - exp(-r * days_hol)
AR = c(traveller_AR, low_AR, high_AR,  endemic_AR)

# add trial populations

all_risk$trial_pop =  c(2736, 346) # from table 1 https://www.sciencedirect.com/science/article/pii/S0140673623006414?via%3Dihub#sec1
all_risk$trial_deaths =  c(0, 0)
all_risk$trial_sae =  c(1, 1)

# Extract beta parameters: default prior is beta(1,1)

# main 
di = extract_beta_params(x = all_risk$deaths, n = all_risk$incidence)
dv = extract_beta_params(x = all_risk$v_deaths, n = all_risk$doses)
ci = extract_beta_params(x = all_risk$cases, n = all_risk$incidence)
cv = extract_beta_params(x = all_risk$SAE, n = all_risk$doses)

# trial data 
dvt = extract_beta_params(x = all_risk$trial_deaths, n = all_risk$trial_pop)
cvt = extract_beta_params(x = all_risk$trial_sae, n = all_risk$trial_pop)

# sensitivity analysis 
ci_half =  extract_beta_params(x = (all_risk$cases/2), n = all_risk$incidence)
dv_3 = extract_beta_params(x = c(0,3), n = all_risk$doses)
dv_0 = extract_beta_params(x = c(0,0), n = all_risk$doses)

# calculate posterior probability benefit 
all_data = simulate_vaccine_benefit(VE, AR, all_risk, di, dv, ci, cv, n, pop) 
trial_data = simulate_vaccine_benefit(VE, AR, all_risk, di, dvt, ci, cvt, n, pop) 
all_data_VE50 = simulate_vaccine_benefit(0.5, AR, all_risk, di, dv, ci, cv, n, pop) 
trial_data_VE50 = simulate_vaccine_benefit(0.5, AR, all_risk, di, dvt, ci, cvt, n, pop) 
all_data_ci_half = simulate_vaccine_benefit(VE, AR, all_risk, di, dv, ci_half, cv, n, pop) 
trial_data_ci_half = simulate_vaccine_benefit(VE, AR, all_risk, di, dvt, ci_half, cvt, n, pop) 
all_data_dv_3 = simulate_vaccine_benefit(VE, AR, all_risk, di, dv_3, ci, cv, n, pop) 
all_data_dv_0 = simulate_vaccine_benefit(VE, AR, all_risk, di, dv_0, ci, cv, n, pop) 


# main figure 
p1 = all_data %>% 
  bind_rows(trial_data, .id = "data") %>% 
  select(data, ar_scenario, age_group, d_benefit, c_benefit) %>% 
  mutate(data = factor(data, labels = c("All data", "Trial data")),
         ar = factor(ar_scenario, levels = c(2,3,4,1),
                     labels = c("Small outbreak",
                                "Large outbreak",
                                "Endemic",
                                "Traveller"))) %>% 
  pivot_longer(c(d_benefit, c_benefit)) %>% 
  mutate(outcome = factor(name, levels = c("c_benefit", "d_benefit"),
                      labels = c("SAE/clinically-attended chikungunya", "Deaths"))) %>% 
  ggplot(aes(x = value, y = age_group)) +
  geom_point(aes(color = ar, shape = data), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", x = "Probability vaccine benefit\noutweighs risk") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  facet_wrap(~outcome)

# VE 0.5 

p2 = all_data_VE50 %>% 
  bind_rows(trial_data_VE50, .id = "data") %>% 
  select(data, ar_scenario, age_group, d_benefit, c_benefit) %>% 
  mutate(data = factor(data, labels = c("All data", "Trial data")),
         ar = factor(ar_scenario, levels = c(2,3,4,1),
                     labels = c("Small outbreak",
                                "Large outbreak",
                                "Endemic",
                                "Traveller"))) %>% 
  pivot_longer(c(d_benefit, c_benefit)) %>% 
  mutate(outcome = factor(name, levels = c("c_benefit", "d_benefit"),
                          labels = c("SAE/clinically-attended chikungunya", "Deaths"))) %>% 
  ggplot(aes(x = value, y = age_group)) +
  geom_point(aes(color = ar, shape = data), size = 2) +
  xlim(0,1)  +
  my_theme + theme(legend.position = "left") +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", x = "Probability vaccine benefit\noutweighs risk") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  facet_wrap(~outcome)

# p case given inf is 2 fold lower  

p3 = all_data_ci_half %>% 
  bind_rows(trial_data_ci_half, .id = "data") %>% 
  select(data, ar_scenario, age_group, c_benefit) %>% 
  mutate(data = factor(data, labels = c("All data", "Trial data")),
         ar = factor(ar_scenario, levels = c(2,3,4,1),
                     labels = c("Small outbreak",
                                "Large outbreak",
                                "Endemic",
                                "Traveller"))) %>% 
  ggplot(aes(x = c_benefit, y = age_group)) +
  geom_point(aes(color = ar, shape = data), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", x = "Probability vaccine benefit\noutweighs risk of SAE") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) 


# 3 vaccine deaths not 1 

p4 = all_data_dv_0 %>% 
  bind_rows(all_data, all_data_dv_3, .id = "data") %>% 
  select(data, ar_scenario, age_group, d_benefit) %>% 
  mutate(data = factor(data, labels = c("0 deaths 65+", "1 death 65+\n(baseline)", "3 death 65+")),
         ar = factor(ar_scenario, levels = c(2,3,4,1),
                     labels = c("Small outbreak",
                                "Large outbreak",
                                "Endemic",
                                "Traveller"))) %>% 
  filter(age_group == "65+") %>% 
  ggplot(aes(x = d_benefit, y = ar)) +
  geom_point(aes(color = data), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5, 8)) +
  labs(y = " ", x = "Probability vaccine benefit\noutweighs risk of death") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) 


# save all plots 

ggsave(
  p1,
  filename = "output/Fig4.pdf",
  units = "cm",
  height = 6,
  width = 14
)

ggsave(
  p2,
  filename = "output/prob_vac_benefit_VE50_new.jpg",
  units = "cm",
  height = 6,
  width = 14
)

ggsave(
  p3,
  filename = "output/prob_vac_benefit_ci_half_new.jpg",
  units = "cm",
  height = 8,
  width = 14
)

ggsave(
  p4,
  filename = "output/prob_vac_benefit_death_sa_new.jpg",
  units = "cm",
  height = 8,
  width = 14
)


# sample size to conclude benefit of non-fatal SAE

cv_n = extract_beta_params(x = c(4,1), n = c(4000,1000))

trial_data = simulate_vaccine_benefit(VE, AR, all_risk, di, dvt, ci, cv_n, n, pop) 

trial_data %>%  select(c_benefit, age_group, ar_scenario)



# add prob plot (from file 7)

FigS7 =  plot_grid(p2, VE_ar,
                   rel_heights = c(1,2),
                   labels = c("a", "b"), ncol = 1)

ggsave(FigS7, 
       filename = "output/SupFig6.jpg",
       units = "cm",
       height = 18,
       width = 14)
