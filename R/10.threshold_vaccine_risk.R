# Script to calculate the maximum acceptable vaccine risk for different diseases, attack rates, vaccine efficacies. 
# File outputs Figure 3.

# Set up 
library(tidyverse)

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

# Import chik data 
chik_risk = read.csv("data/all_risk.csv")[,-1]

chik_risk = chik_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) %>% 
  arrange(age_group) 

# Define parameter ranges 

ar = seq(0,1,0.01) # attack rate
ve = c(0.5, 0.95)
p_symp = seq(0,0.7,0.01)

df = expand.grid(ar, ve, p_symp)
colnames(df) = c("ar", "ve", "p_symp")

# Calculate acceptable vaccine risk for all ar, ve, p_symp
df_v_risk = df %>%
  mutate(v_risk = p_symp * ar * ve) %>% 
  mutate(ar = ar * 100) %>% 
  mutate(ve = paste0("VE = ", ve * 100, "%")) 

# Calculate pathogen specific vaccine risk 


# chik 


chik_risk = chik_risk %>% 
  rename(p_symp = p_case) %>% 
  mutate(ar = 0.3,
         ve_high = 0.95,
         ve_low = 0.5, 
         v_risk_low = p_symp * ve_low * ar,
         v_risk_high = p_symp * ve_high * ar)  %>% 
  pivot_longer(cols = v_risk_low: v_risk_high, values_to = "v_risk") %>% 
  mutate(ve = ifelse(name == "v_risk_low", ve_low, ve_high)) %>%  
  mutate(ve = paste0("VE = ", ve * 100, "%"),
         ar = ar * 100) %>% 
  mutate(disease = "CHIKV")


# covid 19 

daily_ar_2020 = 0.000714
daily_ar_2022 = 0.005713

ar_2020 = (1-(1-daily_ar_2020)^90) # 90 days
ar_2022 = (1-(1-daily_ar_2022)^90) # 90 days

covid_risk = data.frame(
  age = c("25-44", "75+", "25-44", "75+", "25-44", "75+", "25-44", "75+"),
  time = c("2020", "2020", "2022", "2022", "2020", "2020", "2022", "2022"),
  ve = c(0.95, 0.95, 0.95, 0.95, 0.5, 0.5, 0.5, 0.5),
  p_symp = c(0.008, 0.31, 0.002, 0.05, 0.008, 0.31, 0.002, 0.05),
  ar = c(ar_2020, ar_2020, ar_2022, ar_2022, ar_2020, ar_2020, ar_2022, ar_2022)
)

covid_risk = covid_risk %>% 
  mutate(v_risk = p_symp * ve * ar) %>% 
  mutate(ve = paste0("VE = ", ve * 100, "%"),
         ar = ar * 100) %>% 
  mutate(disease = ifelse(time == 2020, "SARS-CoV-2 (2020)", "SARS-CoV-2 (2022)"))

# ebola 

ebola_risk = data.frame(
  ar = c(0.07, 0.07), 
  ve = c(0.5, 0.95),
  p_symp = c(0.6, 0.6)
)


ebola_risk = ebola_risk %>% 
  mutate(v_risk = p_symp * ve * ar) %>% 
  mutate(ve = paste0("VE = ", ve * 100, "%"),
         ar = ar * 100) %>% 
  mutate(disease = "EV")

# measles

measles_risk = data.frame(
  p_symp = c(0.2, 0.2, 0.2, 0.2),
  ve = c(0.5, 0.5, 0.95, 0.95),
  ar = c(0.002, 0.9, 0.002, 0.9),
  disease = c("MeV (now)", "MeV (pre-vaccination)")
)

measles_risk = measles_risk %>% 
  mutate(v_risk = p_symp * ve * ar) %>% 
  mutate(ve = paste0("VE = ", ve * 100, "%"),
         ar = ar * 100) 

# Plot all 

general_fig = df_v_risk %>% 
  ggplot(aes(x = ar, y = p_symp, fill = v_risk)) +
  geom_tile() +
  geom_point(data = covid_risk,
             aes(x = ar, y = p_symp, shape = age, color = disease), size = 3) +
  geom_point(data = chik_risk,
             aes(x = ar, y = p_symp, shape = age_group, color = disease), size = 3) +
  geom_point(data = ebola_risk,
             aes(x = ar, y = p_symp, color = disease), size = 3) +
  geom_point(data = measles_risk,
             aes(x = ar, y = p_symp, color = disease), size = 3) +
  facet_wrap(~ve) +
  my_theme +
  labs(x = "Attack rate (%)", 
       y = "Probability of severe outcome\nfollowing infection", 
       fill = "Risk threshold for\nconcluding vaccine benefit",
       shape = "Age",
       color = "Disease") +
  scale_fill_gradientn(
    colours = c("white", "#E0F6FF", "#FFA500", "red"),
    values = scales::rescale(c(0, 0.05, 0.3, 1))) +
  scale_shape_manual(values = c(
    "25-44"  = 0, 
    "18-64"  = 2,  
    "65+"    = 17,
    "75+"    = 15 
    
  )) +
  scale_color_manual(values = c(
    "CHIKV" = "#9966CC", 
    "EV"       = "#74C476",  
    "SARS-CoV-2 (2020)" = "#CC3399",
    "SARS-CoV-2 (2022)" = "#FFADD6",
    "MeV (now)" = "#00B4D8",
    "MeV (pre-vaccination)" = "#000066"
  )) +
  guides(
    fill  = guide_colorbar(order = 1),  
    color = guide_legend(order = 2),  
    shape = guide_legend(order = 3)    
  )

general_fig

# Save
ggsave(general_fig, 
       filename = "output/general_disease_figure.jpg",
       units = "cm",
       height = 12, width = 20)

  
# Comparing covid est for text

covid_risk %>% 
  arrange(v_risk) %>% 
  select(-ar, -p_symp, -disease) %>%  
  mutate(v_risk=v_risk*100) %>% 
  pivot_wider(names_from = time, values_from = v_risk) 
