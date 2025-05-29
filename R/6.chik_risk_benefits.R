
# Calculate IXCHIQ risks and benefits


library(Hmisc)
library(truncnorm)
library(tidyverse)
library(cowplot)
library(patchwork)

dir.create("output", showWarnings = FALSE)

my_theme = theme_void() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title.y = element_blank(),
    legend.position = "none",
    strip.text = element_blank()
  )

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


# Initialize storage for results

summary_results <- tibble()

# Nested loops over i, j, and k
for (i in 1:length(VE)) {
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
      d_threshold = pdi * VE[i] * AR[k]
      c_threshold = pci * VE[i] * AR[k]
    
      # Calculate proportion averted 
      d_av = (d_threshold - pdv) * pop
      c_av = (c_threshold - pcv) * pop
      
      # Calculate cases with and without vaccination 
      ndi = pdi * AR[k] * pop 
      nci = pci * AR[k] * pop 
      ncv = nci  * (1-VE[i])  
      ndv = ndi  * (1-VE[i]) 
      tdv = (pdi * AR[k] * (1-VE[i]) + pdv) * pop
      tcv = (pci * AR[k] * (1-VE[i]) + pcv) * pop
      
      VE_value = VE[i]
      AR_value = AR[k]
      age_value = all_risk$age_group[j]
   
      # Store summary results
      summary_results <- bind_rows(summary_results, 
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
                                     high_c_av = quantile(c_av, probs = c(0.975)),
                                     mean_pd_i = mean(pdi* pop),
                                     mean_pd_v = mean(pdv* pop),
                                     mean_pc_i = mean(pci* pop),
                                     mean_pc_v = mean(pcv* pop),
                                     low_pd_i  = quantile(pdi* pop, probs = c(0.025)),
                                     low_pd_v  = quantile(pdv* pop, probs = c(0.025)),
                                     high_pd_i = quantile(pdi* pop, probs = c(0.975)),
                                     high_pd_v = quantile(pdv* pop, probs = c(0.975)),
                                     low_pc_i  = quantile(pci* pop, probs = c(0.025)),
                                     low_pc_v  = quantile(pcv* pop, probs = c(0.025)),
                                     high_pc_i = quantile(pci* pop, probs = c(0.975)),
                                     high_pc_v = quantile(pcv* pop, probs = c(0.975)),
                                     mean_nd_i = mean(ndi),
                                     mean_nd_v = mean(ndv),
                                     mean_nc_i = mean(nci),
                                     mean_nc_v = mean(ncv),
                                     low_nd_i  = quantile(ndi, probs = c(0.025)),
                                     low_nd_v  = quantile(ndv, probs = c(0.025)),
                                     high_nd_i = quantile(ndi, probs = c(0.975)),
                                     high_nd_v = quantile(ndv, probs = c(0.975)),
                                     low_nc_i  = quantile(nci, probs = c(0.025)),
                                     low_nc_v  = quantile(ncv, probs = c(0.025)),
                                     high_nc_i = quantile(nci, probs = c(0.975)),
                                     high_nc_v = quantile(ncv, probs = c(0.975)),
                                     mean_td_v = mean(tdv),
                                     mean_tc_v = mean(tcv),
                                     low_td_v  = quantile(tdv, probs = c(0.025)),
                                     low_tc_v  = quantile(tcv, probs = c(0.025)),
                                     high_td_v = quantile(tdv, probs = c(0.975)),
                                     high_tc_v = quantile(tcv, probs = c(0.975)),
                                   ))
          }
  }
}


summary_results = summary_results %>%  
  mutate(scenario = factor(ar_scenario, levels = c(2,3,4,1),
                           labels = c("Small\noutbreak",
                                      "Large\noutbreak",
                                      "Endemic",
                                      "Traveller")))


cases_averted2 = summary_results %>% 
  select(scenario,age_group,  mean_d_av:high_c_av) %>% 
  pivot_longer(-c(age_group, scenario)) %>% 
  separate(name, into = c("stat", "name", NA)) %>% 
  mutate(name = factor(name, levels = c("c", "d"), 
                       labels = c("Severe cases", "Deaths"))) %>% 
  pivot_wider(names_from = stat) %>% 
  mutate(color = ifelse(mean <0, "0", "1")) 
  
p1 = cases_averted2 %>%
  filter(name == "Severe cases") %>%
  ggplot(aes(x = mean, y = age_group)) +
  geom_bar(stat = "identity", aes(fill = color)) +
  geom_errorbar(aes(xmin = low, xmax = high), width =0.25) +
  facet_wrap( ~ scenario , ncol = 1) +
  my_theme +
  scale_fill_manual(values = c("#1B3B6F")) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-400, 400))


p2 = cases_averted2 %>% 
  filter(name != "Severe cases") %>%  
  ggplot(aes(x = mean, y = age_group)) +
  geom_bar(stat = "identity", aes(fill = color)) +
  geom_errorbar(aes(xmin = low, xmax = high),  width =0.25) +
  facet_wrap(~ scenario ,  ncol = 1) +
  my_theme +
  scale_fill_manual(values= c( "#FF6F61", "#1B3B6F")) +
  geom_vline(xintercept = 0, linetype = 2)+
  scale_x_continuous(limits = c(-5,5))

out1 = plot_grid(
  p1, p2, ncol = 1
)

out1

ggsave(
  out1,
  filename = "output/cases_averted_CrI.pdf",
  units = "cm",
  height = 16,
  width = 10,
)

# number of cases --------------------------------------------------------------

cases = summary_results %>% 
  select(scenario,age_group,  mean_nd_i:high_nc_v, 
         mean_pd_v, low_pd_v, high_pd_v,
         mean_pc_v, low_pc_v, high_pc_v,
         mean_tc_v, mean_td_v,
         low_tc_v, high_tc_v, 
         low_td_v, high_td_v) %>% 
  pivot_longer(-c(age_group, scenario)) %>% 
  separate(name, into = c("stat", "name", "cause")) %>% 
  mutate(name = factor(name, levels = c("nc", "nd", "pc", "pd", "tc", "td"), 
                       labels = c("Severe cases", "Deaths", "V Severe cases", "V Deaths", "TV Severe cases", "TV Deaths"))) %>% 
  mutate(cause = factor(cause, levels = c("i", "v"), labels = c("Infection", "Vaccine"))) %>% 
  mutate(value = ifelse(cause == "Infection", -value, value)) %>% 
  pivot_wider(names_from = stat, values_from = value) 

# Left plot
p_left_death = cases %>% 
  filter(name == "Deaths", cause == "Infection") %>% 
  ggplot(aes(y = mean, x = age_group, fill = cause)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low, ymax = high), width=0.25) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Infection" = "#A0BFD8")) +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme +
  scale_y_continuous(labels = abs, limits = c(-5,0)) 


# Right plot

# Filter for mean deaths attributal to vaccine or infection 
death_mean_data = cases %>% 
  filter(name %in% c("V Deaths", "Deaths"), cause != "Infection") %>% 
  select(-low, -high)

# Filter for error bars (from total deaths only (vaccine + infection))
death_error_data = cases %>% 
  filter(name == "TV Deaths", cause != "Infection") %>% 
  select(age_group, scenario, low, high)

# Merge error bars into mean data
death_plot_data = death_mean_data %>% 
  left_join(death_error_data, by = c("age_group", "scenario"))


p_right_death = death_plot_data %>% 
  ggplot(aes(y = mean, x = age_group, fill = name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low, ymax = high), width=0.25) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("V Deaths" = "#A3B18A",
                               "Deaths" = "#4B4B4B")) +
  my_theme +
  facet_wrap(~ scenario, ncol = 1) +
  labs(y = "Deaths (per 10,000)") +
  scale_y_continuous(labels = abs, limits = c(0,6))

# Left plot
p_left_dis = cases %>% 
  filter(name == "Severe cases", cause == "Infection") %>%  
  ggplot(aes(y = mean, x = age_group, fill = cause)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low, ymax = high), width=0.25) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Infection" = "#A0BFD8")) +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme +
  scale_y_continuous(labels = abs, limits = c(-400,0)) 


# Right plot

# Filter for mean cases attributal to vaccine or infection 
cases_mean_data = cases %>% 
  filter(name %in% c("V Severe cases", "Severe cases"), cause != "Infection") %>% 
  select(-low, -high)

# Filter for error bars (from total cases only (vaccine + infection))
cases_error_data = cases %>% 
  filter(name == "TV Severe cases", cause != "Infection") %>% 
  select(age_group, scenario, low, high)

# Merge error bars into mean data
case_plot_data = cases_mean_data %>% 
  left_join(cases_error_data, by = c("age_group", "scenario"))

p_right_dis = case_plot_data %>% 
  ggplot(aes(y = mean, x = age_group, fill = name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low, ymax = high), width=0.25) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("V Severe cases" = "#A3B18A",
                               "Severe cases" = "#4B4B4B")) +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme +
  scale_y_continuous(labels = abs, limits = c(0,400))


# combine plots

# Combine the plots
combined_plot_death = plot_grid(p_left_death, p_right_death)

combined_plot_disease = plot_grid(p_left_dis, p_right_dis)

ggsave(combined_plot_disease/combined_plot_death,
       filename = "output/all_with_without_vaccination_CrI.pdf",
       units = "cm",
       height = 16,
       width = 10)

