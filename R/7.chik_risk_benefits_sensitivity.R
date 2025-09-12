# Script to run sensitivity analyses on the impact of VE and AR on expected number of medically 
# attended chikungunya cases / SAEs and deaths averted by IXCHIQ vaccination. 
# File outputs Figures S4 and S5. 

library(Hmisc)
library(truncnorm)
library(tidyverse)
library(cowplot)
library(patchwork)

dir.create("output", showWarnings = FALSE)


my_theme = theme_bw() +
  theme(
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )

all_risk = read.csv("data/all_risk.csv")[,-1]

all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) %>% 
  mutate(age_group = fct_rev(age_group))

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


################################################################################
########################## FOI sensitivity analysis ############################
################################################################################

summary_results_AR = tibble()

AR_SA = seq(0,0.5,0.01)

# Nested loops over i, j, and k
for (i in 1:length(VE)) {
  for (j in 1:nrow(all_risk)) {
    for (k in 1:length(AR_SA)) {
    
      
      # Generate random samples
      ndi = rbinom(n, size = round(all_risk$incidence)[j], prob = all_risk$ifr[j])
      ndv = rbinom(n, size = round(all_risk$doses)[j], prob = all_risk$p_deaths[j])
      nci = rbinom(n, size = round(all_risk$incidence)[j], prob = all_risk$p_case[j])
      ncv = rbinom(n, size = round(all_risk$doses)[j], prob = all_risk$p_sae[j])
      
      pdi = ndi / round(all_risk$incidence)[j]
      pdv = ndv / round(all_risk$doses)[j]
      pci = nci / round(all_risk$incidence)[j]
      pcv = ncv / round(all_risk$doses)[j]
      
      # Calculate threshold and prob benefit
      d_threshold = pdi * VE[i] * AR_SA[k]
      c_threshold = pci * VE[i] * AR_SA[k]
    
      # Calculate proportion averted 
      d_av = (d_threshold - pdv) * pop
      c_av = (c_threshold - pcv) * pop
    
      VE_value = VE[i]
      AR_value = AR_SA[k]
      age_value = all_risk$age_group[j]
   
      # Store summary results
      summary_results_AR = bind_rows(summary_results_AR, 
                                   tibble(
                                     ar_scenario = k, 
                                     age_group = age_value,
                                     VE_value = VE_value,
                                     AR_value = AR_value,
                                     mean_d_av = mean(d_av),
                                     mean_c_av = mean(c_av),
                                     low_d_av  = quantile(d_av, probs = c(0.025)),
                                     low_c_av  = quantile(c_av, probs = c(0.025)),
                                     high_d_av = quantile(d_av, probs = c(0.975)),
                                     high_c_av = quantile(c_av, probs = c(0.975))
                                   ))
          }
  }
}

ar_df = data.frame(
  ar_values = AR*100,
  # Create two rows for each AR value - one for each facet
  name = rep(factor(c("Severe cases", "Deaths")), each = length(AR))
)


cases_averted_AR = summary_results_AR %>% 
  select(AR_value, age_group,  mean_d_av:high_c_av) %>% 
  pivot_longer(-c(age_group, AR_value)) %>% 
  separate(name, into = c("stat", "name", NA)) %>% 
  mutate(name = factor(name, levels = c("c", "d"), 
                       labels = c("Severe cases", "Deaths"))) %>% 
  pivot_wider(names_from = stat) 
  
p1 = cases_averted_AR %>%
  filter(name == "Deaths") %>% 
  ggplot(aes(x = AR_value*100, y = mean)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = age_group), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  mapply(function(ar_val) {
    list(
      annotate("segment", 
               x = ar_val*100, xend = ar_val*100,
               y = -Inf, yend = -1.5,
               arrow = arrow(length = unit(0.1, "cm")))
    )
  }, AR) +
  labs(x = "Attack rate (%)",
       y = "Deaths averted\n(per 10,000 vaccinated)") +
    my_theme +
  theme(legend.position = "none") +
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  scale_fill_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  scale_y_continuous(breaks = c(-2, 0, 2, 4, 6))



p2 = cases_averted_AR %>%
  filter(name != "Deaths") %>% 
  ggplot(aes(x = AR_value*100, y = mean)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = age_group), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  mapply(function(ar_val) {
    list(
      annotate("segment", 
               x = ar_val*100, xend = ar_val*100,
               y = -Inf, yend = -40,
               arrow = arrow(length = unit(0.1, "cm")))
    )
  }, AR) +
  labs(x = " ",
       y = "Clinically-attented chikungunya\ncases / SAEs averted\n(per 10,000 vaccinated)") +
  my_theme +
  theme(legend.position = c(0.07,0.76),
        legend.title = element_blank()) +
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  scale_fill_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  scale_y_continuous(breaks = c(-40, 0, 200, 400, 600))

SA_ar = p2/p1

SA_ar

ggsave(SA_ar, 
       filename = "output/AR_SA_plot_new.jpg",
       units = "cm",
       height = 10,
       width = 14)

################################################################################
########################## VE sensitivity analysis #############################
################################################################################

summary_results_VE= tibble()

VE_sa = seq(0,1,0.01)

# Nested loops over i, j, and k
for (i in 1:length(VE_sa)) {
  for (j in 1:nrow(all_risk)) {
    for (k in 1:length(AR)) {
      
      
      # Generate random samples
      ndi = rbinom(n, size = round(all_risk$incidence)[j], prob = all_risk$ifr[j])
      ndv = rbinom(n, size = round(all_risk$doses)[j], prob = all_risk$p_deaths[j])
      nci = rbinom(n, size = round(all_risk$incidence)[j], prob = all_risk$p_case[j])
      ncv = rbinom(n, size = round(all_risk$doses)[j], prob = all_risk$p_sae[j])
      
      pdi = ndi / round(all_risk$incidence)[j]
      pdv = ndv / round(all_risk$doses)[j]
      pci = nci / round(all_risk$incidence)[j]
      pcv = ncv / round(all_risk$doses)[j]
      
      # Calculate threshold and prob benefit
      d_threshold = pdi * VE_sa[i] * AR[k]
      c_threshold = pci * VE_sa[i] * AR[k]
      
      # Calculate proportion averted 
      d_av = (d_threshold - pdv) * pop
      c_av = (c_threshold - pcv) * pop
      
      VE_value = VE_sa[i]
      AR_value = AR[k]
      age_value = all_risk$age_group[j]
      
      # Store summary results
      summary_results_VE = bind_rows(summary_results_VE, 
                                     tibble(
                                       ar_scenario = k, 
                                       age_group = age_value,
                                       VE_value = VE_value,
                                       AR_value = AR_value,
                                       mean_d_av = mean(d_av),
                                       mean_c_av = mean(c_av),
                                       low_d_av  = quantile(d_av, probs = c(0.025)),
                                       low_c_av  = quantile(c_av, probs = c(0.025)),
                                       high_d_av = quantile(d_av, probs = c(0.975)),
                                       high_c_av = quantile(c_av, probs = c(0.975))
                                     ))
    }
  }
}

cases_averted_VE = summary_results_VE %>% 
  select(VE_value, ar_scenario, age_group,  mean_d_av:high_c_av) %>% 
  pivot_longer(-c(age_group, VE_value, ar_scenario, VE_value)) %>% 
  separate(name, into = c("stat", "name", NA)) %>% 
  mutate(name = factor(name, levels = c("c", "d"), 
                       labels = c("Severe cases", "Deaths"))) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(ar_scenario = factor(ar_scenario, levels = c(2,3,4,1),
                              labels = c("Small outbreak",
                                         "Large outbreak",
                                         "Endemic",
                                         "Traveller"))) 

p3 = cases_averted_VE %>%
  filter(name == "Deaths") %>% 
  ggplot(aes(x = VE_value*100, y = mean)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = age_group), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 95, linetype = 2, colour = "black") +
  labs(x = "Vaccine efficacy (%)", 
       y = "Deaths averted\n(per 10,000 vaccinated)") +
  my_theme +
  theme(legend.position ="none") +
  scale_x_continuous(limits  = c(0,100))+
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  scale_fill_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  facet_wrap(~ar_scenario, ncol=4)

p4 = cases_averted_VE %>%
  filter(name != "Deaths") %>% 
  ggplot(aes(x = VE_value*100, y = mean)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = age_group), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 95, linetype = 2, colour = "black") +
  labs(x = "", 
       y = "Clinically-attented chikungunya\ncases / SAEs averted\n(per 10,000 vaccinated)") +
  my_theme +
  theme(legend.position = c(0.07,0.76),
        legend.title = element_blank()) +
  scale_x_continuous(limits  = c(0,100))+
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  scale_fill_manual(values = c(
    "#FD8D3C", 
    "#C51B8A"
  )) +
  facet_wrap(~ar_scenario, ncol=4) +
  scale_y_continuous(breaks = c(-30, 0, 200, 400))



VE_ar = p4/p3

VE_ar

ggsave(VE_ar, 
       filename = "output/VE_SA_plot_new.jpg",
       units = "cm",
       height = 10,
       width = 14)

