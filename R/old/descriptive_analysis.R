
library(tidyverse)
library(patchwork)
# Calculate risks and plot 

# set common theme 
my_theme = theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

# --- import data

inf_data = read.csv("data/processed_data.csv")[,-1]
demog = read.csv("data/reunion_population_age_groups_2021.csv")

#  ---- define other data 
doses_data = tibble(
  age_group = c("<65", "65+"),
  doses = c(21282, 16236)
)

doses_data_2 = tibble(
  age_group = c("18-64", "65+", "18-64", "65+", "18-64", "65+", "18-64", "65+"),
  location = c("France Mainland\n& Overseas", "France Mainland\n& Overseas", "La Réunion",  "La Réunion",
               "US", "US", "Canada", "Canada"),
  doses = c(1206, 4823, 824, 5594, 14675, 4675, 4577, 1144)
)

write.csv(doses_data_2, "data/doses_country.csv")


# SAE by age groups 
vaccine_risk = tibble(
  age_group = c("18-64", "65-69", "70-79", "80+"),
  SAE = c(3, 3, 3, 6),
  v_deaths = c(0, 0, 1, 1)
)

# --- get doses by plus 65 categories

plus_65_demo = demog %>%  
  slice(14:21) %>% 
  mutate(age_group = c("65-69","70-79", "70-79", "80+","80+","80+","80+","80+")) %>% 
  group_by(age_group) %>% 
  summarise(total = sum(Total)) %>% 
  mutate(prop = total / sum(total))

write.csv(plus_65_demo, "data/plus_65_demo.csv")

total_65_doses = doses_data %>% 
  filter(age_group == "65+") %>% 
  pull(doses)

plus_65_doses = plus_65_demo %>% 
  mutate(total_doses = total_65_doses,
         doses = total_doses * prop) 

vaccine_risk = vaccine_risk %>% 
  mutate(doses = c(doses_data$doses[1], plus_65_doses$doses)) %>% 
  mutate(p_sae = SAE/doses) %>% 
  mutate(p_deaths = v_deaths/doses)
  
# combine inf and vaccine risk 
all_risk = inf_data %>%  
  left_join(vaccine_risk) %>% 
  mutate(age_group = factor(age_group, levels = c("18-64", "65-69", "70-79", "80+")))

# Map consistent values to the measure column

death_risk_data = all_risk %>%
  pivot_longer(cols = c(ifr, p_deaths), names_to = "measure", values_to = "value") %>%
  mutate(measure = recode(measure,
                          "ifr" = "Infection",
                          "p_deaths" = "Vaccine"))

disease_risk_data = all_risk %>%
  pivot_longer(cols = c(p_case, p_sae), names_to = "measure", values_to = "value") %>%
  mutate(measure = recode(measure,
                          "p_case" = "Infection",
                          "p_sae" = "Vaccine"))

# Combine data for a unified legend
common_palette = c("Vaccine" = "#FDB863",
                   "Infection" = "#49006A")

# Death risk bar plot
death_risk = death_risk_data %>% 
  ggplot(aes(x = age_group, y = value, fill = measure)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 0.004)) +
  scale_fill_manual(values =common_palette) +
  labs(y = "Probability of\ndeath",
       x = "Age group") +
  my_theme +
  theme(legend.position = c(0.02, 0.98),  # top-left
        legend.justification = c(0, 1))


# disease risk bar plot with legend
disease_risk = disease_risk_data %>% 
  ggplot(aes(x = age_group, y = value, fill = measure)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_fill_manual(values =common_palette) +
  labs(y = "Probability of\nclinically-attented\nchikungunya / SAEs",
       x = "Age group") +
  my_theme +
  theme(legend.position = "none")

risk_plots = disease_risk / death_risk +
  plot_layout(guides = "collect")

risk_plots

ggsave(risk_plots, 
       filename = "output/Fig1b.jpeg",
       units = "cm",
       height = 12, width = 15)


write.csv(all_risk, "data/all_risk.csv")



# First, extract the individual plots from each patchwork


# Plot 1: Doses by age group
doses_plot = ggplot(doses_data_2, aes(x = age_group, y = doses)) +
  geom_bar(stat = "identity", aes(fill = location), width = 0.7) +
  labs(
    y = "Vaccine doses",
    x = "Age group"
  ) +
  ylim(0,25000) +
  my_theme +
  theme(legend.position = "top") +
  scale_fill_manual(values = c(
    "#FD8D3C", 
    "#F768A1",
    "#C51B8A",
    "#7A0177"  
  )) +
  theme(
    legend.position = c(0.95, 0.95),  # Position in the top-right (x=0.95, y=0.95)
    legend.justification = c(1, 1),    # Anchor at the top-right
    legend.direction = "horizontal",   # Make legend horizontal
    legend.background = element_rect(fill = "white", color = NA, linewidth = 0),
    legend.key.size = unit(0.4, "cm"),
    legend.key.width = unit(0.4, "cm")
  )


# Plot 2: SAE by age group
SAE_plot = ggplot(vaccine_risk, aes(x = age_group, y = SAE)) +
  geom_bar(stat = "identity", fill = "#FDB863", width = 0.7) +
  labs(
    y = "Vaccine assocaited\nSAEs",
    x = "Age group"
  ) +
  ylim(0,6) +
  my_theme


# Plot 3: Deaths by age group
deaths_plot = ggplot(vaccine_risk, aes(x = age_group, y = v_deaths)) +
  geom_bar(stat = "identity", fill = "#FDB863", width = 0.7) +
  labs(
    y = "Vaccine associated\ndeaths",
    x = "Age group"
  ) +
  ylim(0,6) +
  my_theme


# 
# p1 = (SAE_plot | deaths_plot | doses_plot)+
#   plot_layout(guides = "collect")
# 
# p2  = (diSAEse_risk | death_risk) +
#   plot_layout(guides = "collect") +
#   theme(legend.position = c(0.5,0.9)) 
# 
# p3 = (p1 / p2)  +
#   plot_annotation(tag_levels = "a") + 
#   plot_layout(axis_titles = "collect")
# 
# p3
# 
# 
# ggsave(p3, 
#        filename = "output/Fig1_V1.jpeg",
#        units = "cm",
#        height = 12, width = 20)
# 


col1 = SAE_plot / disease_risk
col2 = deaths_plot / death_risk

# Final figure
p_final = doses_plot / (col1 | col2) +
  plot_annotation(tag_levels = list(c("a", "b", "d", "c", "e"))) +
  plot_layout(heights = c(1, 2.2))  


p_final

# Save figure
ggsave(p_final, 
       filename = "output/Fig1_V1_updated.jpeg",
       units = "cm",
       height = 14, width = 18)

