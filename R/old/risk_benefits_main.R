
library(tidyverse)
library(cowplot)
library(patchwork)

# Calculate death risk

# --- import data

all_risk = read.csv("data/all_risk.csv")[,-1]

all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) %>% 
  arrange(age_group)

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

population = 10000

results_list = list()

for (i in 1:length(AR)) {


# Create a copy of the original dataframe for this iteration
temp_risk = all_risk


# Calculate risk metrics using the current AR value
temp_risk$risk_death_nv = temp_risk$ifr * AR[i] * population
temp_risk$risk_absdeath_v = temp_risk$p_deaths * population
temp_risk$risk_dis_nv = temp_risk$p_case * AR[i] * population 
temp_risk$risk_absdis_v = temp_risk$p_sae * population

temp_risk$risk_dis_v = (temp_risk$p_case * AR[i]  * (1-VE) + temp_risk$p_sae) * population
temp_risk$risk_death_v = (temp_risk$ifr * AR[i]  * (1-VE) + temp_risk$p_deaths) * population


temp_risk$deaths_averted = temp_risk$risk_death_nv - temp_risk$risk_death_v
temp_risk$cases_averted = temp_risk$risk_dis_nv - temp_risk$risk_dis_v


# Add daily_ar value as a column for reference
temp_risk$daily_ar = AR[i]
temp_risk$scenario = i

# Store the results for this iteration
results_list[[i]] = temp_risk
}

all_results = do.call(rbind, results_list)

all_results_to_plot = all_results %>% 
  select(age_group, risk_death_nv:scenario) %>% 
  mutate(scenario = factor(scenario, levels = c(2,3,4,1),
                           labels = c("Small outbreak",
                                      "Large outbreak",
                                      "Endemic",
                                      "Traveller"))) %>% 
  mutate(risk_chikdeath_v = risk_death_v - risk_absdeath_v,
         risk_chikdis_v = risk_dis_v - risk_absdis_v) %>% 
  pivot_longer(c(risk_death_nv:risk_death_v, c(risk_chikdeath_v, risk_chikdis_v)) , names_to = "outcome", values_to = "N") %>% 
  separate(outcome, into = c(NA, "outcome", "cause"), sep = "_") %>% 
  mutate(N_modified = ifelse(cause == "nv", -N, N)) %>% 
  mutate(cause = factor(cause, levels = c("nv", "v"), labels = c("Infection", "Vaccine")))

# plot pyramid 

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
    strip.text = element_blank(),
  )

# Left plot
p_left_death = all_results_to_plot %>% 
  filter(outcome == "death", cause == "Infection") %>% 
  ggplot(aes(y = N_modified, x = age_group, fill = cause)) +
  geom_bar(stat = "identity") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Infection" = "#A0BFD8")) +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme +
  scale_y_continuous(labels = abs, limits = c(-4,0)) 
  

# Right plot
p_right_death = all_results_to_plot %>% 
  filter(outcome %in% c("absdeath", "chikdeath"), cause != "Infection") %>%  
  ggplot(aes(y = N_modified, x = age_group, fill = outcome)) +
  geom_bar(stat = "identity") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("absdeath" = "#A3B18A",
                               "chikdeath" = "#4B4B4B")) +
  my_theme +
  facet_wrap(~ scenario, ncol = 1) +
  labs(y = "Deaths (per 10,000)") +
  scale_y_continuous(labels = abs, limits = c(0,4))
  
# Left plot
p_left_dis = all_results_to_plot %>% 
  filter(outcome == "dis", cause == "Infection") %>%  
  ggplot(aes(y = N_modified, x = age_group, fill = cause)) +
  geom_bar(stat = "identity") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Infection" = "#A0BFD8")) +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme +
  scale_y_continuous(labels = abs, limits = c(-400,0)) 
  

# Right plot
p_right_dis = all_results_to_plot %>% 
  filter(outcome %in% c("absdis", "chikdis"), cause != "Infection") %>%  
  ggplot(aes(y = N_modified, x = age_group, fill = outcome)) +
  geom_bar(stat = "identity") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("absdis" = "#A3B18A",
                               "chikdis" = "#4B4B4B")) +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme +
  scale_y_continuous(labels = abs, limits = c(0,400)) 
  

# combine plots

# Combine the plots
combined_plot_death = plot_grid(p_left_death, p_right_death)

combined_plot_disease = plot_grid(p_left_dis, p_right_dis)

ggsave(combined_plot_disease/combined_plot_death,
       filename = "output/all_with_without_vaccination.pdf",
       units = "cm",
       height = 16,
       width = 10)


# plot deaths / cases averted --------------------------------------------------
cases_averted = all_results %>% 
  select(age_group, deaths_averted, cases_averted, scenario) %>% 
  mutate(scenario = factor(scenario, levels = c(2,3,4,1),
                           labels = c("Small\noutbreak",
                                      "Large\noutbreak",
                                      "Endemic",
                                      "Traveller"))) %>% 
  pivot_longer(c(deaths_averted, cases_averted)) %>% 
  mutate(name = factor(name, labels = c("Severe cases", "Deaths"))) %>% 
  mutate(color = ifelse(value <0, "0", "1")) 


p1 = cases_averted %>%
  filter(name == "Severe cases") %>%
  ggplot(aes(x = value, y = age_group)) +
  geom_bar(stat = "identity", aes(fill = color)) +
  facet_wrap( ~ scenario , ncol = 1) +
  my_theme +
  scale_fill_manual(values = c("#1B3B6F")) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-400, 400))


  p2 = cases_averted %>% 
    filter(name != "Severe cases") %>%  
    ggplot(aes(x = value, y = age_group)) +
    geom_bar(stat = "identity", aes(fill = color)) +
    facet_wrap(~ scenario ,  ncol = 1) +
    my_theme +
    scale_fill_manual(values= c( "#FF6F61", "#1B3B6F")) +
    geom_vline(xintercept = 0, linetype = 2)+
    scale_x_continuous(limits = c(-2,2))

out1 = plot_grid(
    p1, p2, ncol = 1
  )

  ggsave(
    out1,
    filename = "output/main_figure.pdf",
    units = "cm",
    height = 16,
    width = 10,
  )

