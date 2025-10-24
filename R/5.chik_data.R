# Script to plot the available data on IXCHIQ doses, SAEs, and risk of severe 
# outcomes given chikungunya infection. File outputs Figure S3. 

library(tidyverse)
library(Hmisc)
library(patchwork)

dir.create("output", showWarnings = FALSE)

# set common theme 
my_theme = theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.margin = margin(-3, 0 , -20, -2),
    legend.box.margin = margin(-2, 0 , -2, 0),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

common_palette = c("#A0BFD8", "#FDB863")

# --- import data

all_risk = read.csv("data/all_risk.csv")[,-1]
doses_data_country = read.csv("data/doses_data_us.csv")

# sample 

all = data.frame(
  pdeath_i = binconf(x = all_risk$deaths, n = all_risk$incidence, method = "exact"),
  pcase_i =  binconf(x = all_risk$cases, n = all_risk$incidence, method = "exact"),
  pdeath_v = binconf(x = all_risk$v_deaths, n = all_risk$doses, method = "exact"),
  pcase_v  = binconf(x = all_risk$SAE, n = all_risk$doses, method = "exact")
)

all_long = all %>% 
  mutate(age_group = c("18-64", "65+")) %>% 
  pivot_longer(cols = -age_group) %>% 
  separate(name, into = c("outcome", "cause"), sep = "_") %>% 
  separate(cause, into = c("cause", "stat")) %>% 
  mutate(outcome = factor(outcome, levels = c("pcase", "pdeath"), 
                          labels = c("Medically-attended chikungunya cases\n/ SAEs", "Deaths"))) %>% 
  mutate(cause = factor(cause, levels = c("i", "v"), labels = c("Infection", "Vaccine"))) %>% 
  pivot_wider(names_from = stat)


# Plot 1: Doses by age group
p1 = ggplot(doses_data_country, aes(x = age_group, y = doses)) +
  geom_bar(stat = "identity", aes(fill = location), width = 0.7) +
  labs(
    y = "Vaccine doses",
    x = "Age group"
  ) +
  my_theme +
  theme(legend.position = "top") +
  scale_fill_manual(values = c(
    "#FD8D3C", 
    "#C51B8A",
    "#7A0177"  
  )) +
  theme(
    legend.position = "top",  # Position in the top-right (x=0.95, y=0.95)
    legend.key.size = unit(0.2, "cm")
  )


# Plot 2: SAE by age group
p2 = all_risk %>%
  pivot_longer(cols = c(v_deaths, SAE)) %>%
  mutate(name = factor(name, levels = c("SAE","v_deaths"),
                       labels = c("SAE", "Deaths"))) %>% 
  ggplot(aes(x = age_group, y = value)) +
  geom_bar(stat = "identity", fill = "#FDB863", width = 0.7) +
  labs(
    y = "Vaccine assocaited\nadverse events",
    x = "Age group"
  ) +
  my_theme +
  facet_wrap(~name)

# risk bar plot
p3 = all_long %>% 
  ggplot(aes(x = age_group, y = PointEst, fill = cause)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, group = cause), position = position_dodge(width = 0.8), width  = 0.2) +
  scale_fill_manual(values = common_palette) +
  labs(y = "Probability \nof outcome",
       x = "Age group") +
  my_theme +
  theme(legend.position = c(0.58, 0.92),  # top-left
        legend.justification = c(0, 1)) +
  facet_wrap(~outcome, scale = "free_y")

# Final figure
p_final = (p1 + p2) / p3 + 
  plot_annotation(tag_levels = "a") 

p_final

# Save figure
ggsave(p_final, 
       filename = "output/chik_data.jpg",
       units = "cm",
       height = 10, width = 14)

# get proportions for paper 

all_long %>%  filter(cause == "Vaccine") %>%  mutate(x = paste0(round(PointEst*10000, 1), " (95% CI:", round(Lower*10000,1), " - ", round(Upper*10000,1),")"))


all_long %>%  filter(cause != "Vaccine") %>%  mutate(x = paste0(round(PointEst*10000, 2), " (95% CI:", round(Lower*10000,1), " - ", round(Upper*10000,1),")"))

