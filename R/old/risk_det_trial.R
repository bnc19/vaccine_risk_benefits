### Calculate power to detect adverse events in trial 

library(tidyverse)
library(patchwork)

my_theme = theme_bw() +
  theme(
    legend.text = element_text(size = 7),
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )


# --- import data
all_risk = read.csv("data/all_risk.csv")[,-1]
plus_65_demo = read.csv("data/plus_65_demo.csv")[,-1]

all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65-69", "70-79", "80+"))) %>% 
  arrange(age_group)

# define observed risks 
sae_points = data.frame(
  p_sae = all_risk$p_sae * 100000,
  pop = c(2736, rep(346, 3)), # from table 1 https://www.sciencedirect.com/science/article/pii/S0140673623006414?via%3Dihub#sec1
  age_group = as.factor(all_risk$age_group),
  demo = c(1, plus_65_demo$prop) # break down over 65 to analysis age groups
)

sae_points$pop = sae_points$pop * sae_points$demo


death_points = data.frame(
  p_death = all_risk$p_deaths[c(3, 4)] * 100000,
  pop = sae_points$pop[c(3,4)], 
  age_group = as.factor(all_risk$age_group[c(3, 4)])
)

death_points$age_group = fct_rev(death_points$age_group)
sae_points$age_group = fct_rev(sae_points$age_group)


# Define population sizes
max_n = 35000
population = seq(1, max_n, 100) 

p_death = seq(0, 0.001, 0.00001)
p_sae = seq(0, 0.003, 0.00002)

# Function to calculate power
calculate_power = function(n, p1) {
  power = 1 - pbinom(0, size = n, prob = p1)
  return(power)
}

# Create result matrices
n_deaths = matrix(nrow = length(population), ncol = length(p_death))
n_cases = matrix(nrow = length(population), ncol = length(p_sae))

# Fill n_deaths matrix
for (i in 1:length(p_death)) {
  for (j in 1:length(population)) {
    n_deaths[j, i] = calculate_power(population[j], p_death[i])
  }
}

for (i in 1:length(p_sae)) {
  for (j in 1:length(population)) {
    n_cases[j, i] = calculate_power(population[j], p_sae[i])
  }
}

# format and plot deaths 

df_deaths = as.data.frame(n_deaths)
colnames(df_deaths) = as.character(p_death * 100000)
df_deaths$pop = population

df_deaths_long = df_deaths %>%
  pivot_longer(cols = -pop, names_to = "p_death", values_to = "power") %>%
  mutate(p_death = as.numeric(p_death))


death_plot = df_deaths_long %>% 
  ggplot(aes(x = p_death, y = pop)) + 
  geom_tile(aes(fill = power)) +
  geom_contour(aes(z = power), breaks = 0.8, color = "black", linetype = 2) +
  scale_fill_gradient(
    low = "white", 
    high = "skyblue", 
    limits = c(0, 1)
  ) + 
  my_theme  + 
  geom_point(
    data = death_points,
    aes(x = p_death, y = pop, shape = age_group),
    color = "black",
    size = 2
  ) +
  scale_shape_manual(values = c(3, 4)) +
  labs(x = "Risk of vaccine-associated death (per 100,000)", 
       y = "Trial population", 
       fill = "Power") +
  theme(legend.position =  "none")

# cases ---

df_cases = as.data.frame(n_cases)
colnames(df_cases) = as.character(p_sae * 100000)  # name columns by p_sae in per 100,000
df_cases$pop = population  # add population column

# Convert to long format
df_cases_long = df_cases %>%
  pivot_longer(cols = -pop, names_to = "p_sae", values_to = "power") %>%
  mutate(p_sae = as.numeric(p_sae))  # convert back to numeric
  
case_plot = df_cases_long %>% 
  ggplot(aes(x = p_sae, y = pop)) + 
  geom_tile(aes(fill = power)) +
  geom_contour(aes(z = power), breaks = 0.8, color = "black", linetype = 2) +
  scale_fill_gradient(
    low = "white", 
    high = "skyblue", 
    limits = c(0, 1)
  ) + 
  my_theme  +
  geom_point(
  data = sae_points,
  aes(x = p_sae, y = pop, shape = age_group),
  color = "black",
  size = 2) +
  scale_shape_manual(values = c(16, 17, 3, 4))  +
  labs(x = "Risk vaccine-associated severe disease (per 100,000)", 
       y = "Trial population", 
       fill = "Power",
       shape = "Age group")


# save --- 

out = case_plot / death_plot +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect")

out

ggsave(out,
       filename = "output/risk_det_trial.jpeg",
       units = "cm",
       height = 14, width = 18)

