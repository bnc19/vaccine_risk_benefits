# Script to calculate sample size needed to achieve 80% power to detect at least 
# one death as a function of attack rate and probability of death from infection. 
# File outputs Figure S5. 

library(tidyverse)
library(patchwork)
library(stats)

dir.create("output", showWarnings = FALSE)

my_theme = theme_classic() +
  theme(
    legend.text = element_text(size = 6),
    legend.margin = margin(0, 0, 0, 0),  # Reduced to zero
    legend.box.margin = margin(0, 0, 0, 0),  # Reduced to zero
    legend.spacing.y = unit(0, "cm"),  # Reduced to zero
    legend.spacing.x = unit(0, "cm"),  # Added to control horizontal spacing
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    
  )


# --- import data
all_risk = read.csv("data/all_risk.csv")[,-1]

all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) %>% 
  arrange(age_group)


# Pull out IFR for chik, sars, nipah 

chik = all_risk$ifr[1]
sars = 0.01 
ebola = 0.6

ifr_data = data.frame(
ifr =  c(chik, sars, ebola),
path =  c("CHIKV 18-64 years\n", "SARS-CoV-2", "EV")
)

calculate_power = function(n, p1) {
  power = 1 - pbinom(0, size = n, prob = p1)
  return(power)
}

find_min_n = function(p1, target_power, max_n) {
  for (n in 1:max_n) {
    if (calculate_power(n, p1) >= target_power) {
      return(n)
    }
  }
  stop("No value of n found that achieves the target power within the range.")
}


# Function to calculate sample size for vaccine safety test
calculate_sample_size_safety = function(p_d_I,
                                        lambda_FOI,
                                        VE = 0.95,
                                        power = 0.80,
                                        max_n = 1000) {
  # p0: Maximum acceptable adverse event rate
  p0 = VE * lambda_FOI * p_d_I 
  
  n = find_min_n(p0, target_power = power, max_n = max_n)

  return(n)
}

# Create a grid of values
p_d_I_values = seq(0.00005, 0.61, 0.005)  # Disease severity
lambda_FOI_values = seq(0.01, 0.5, 0.01)  # Infection probability
VE = c(0.95)

# Calculate for different FOI and criterion combinations
results = expand.grid(
  p_d_I = p_d_I_values, 
  lambda_FOI = lambda_FOI_values,
  VE = VE
)

results$n = mapply(calculate_sample_size_safety,
                   p_d_I = results$p_d_I,
                   lambda_FOI = results$lambda_FOI,
                   VE = results$VE,
                   max_n = 4000000)


max(results$n)
min(results$n)

VE_high = results %>%
  filter(VE == 0.95) %>% 
  mutate(n = ifelse(n>20000, 20000, n)) %>% 
  ggplot(aes(x = p_d_I, y = lambda_FOI*100, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    limits = c(1, 20000),  # focus scale up to 5,000
    breaks = c(2, 20, 200, 2000, 20000),
    labels = c(2, 20, 200, 2000, ">20,000")) +
  geom_segment(
    data = ifr_data,
    aes(x = ifr , xend = ifr, y = 0, yend = 2),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  my_theme +
  labs(
    x = "Probability of death following infection",
    y = "Attack rate (%)",
    fill = "Sample size"
  ) +
  theme(
    legend.position = "top",
    legend.box.spacing = unit(0, "cm"),  # Space between legend box and plot
    plot.margin = margin(2, 2, 0, 2)  # Reduce bottom margin of main plot
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(6, "cm"),
      barheight = unit(0.1, "cm"),
      # Add these to control legend bar spacing
      frame.colour = NA,  # Remove frame
      ticks.colour = NA  # Remove ticks
    )
  ) 

# Calculate for different pathogens
results_2 = expand.grid(
  p_d_I = ifr_data$ifr, 
  lambda_FOI = seq(0.01,0.5, 0.01),
  VE = VE
)

results_2$n = mapply(calculate_sample_size_safety,
                   p_d_I = results_2$p_d_I,
                   lambda_FOI = results_2$lambda_FOI,
                   VE = results_2$VE,
                   max_n = 4000000)

risk_2 =  paste0(ifr_data$path , " (", round(ifr_data$ifr*100,3), "%)")

hline_data = results_2 %>%  
  slice(13:15) %>% # check this is is AR = 5% 
  mutate(risk = factor(p_d_I, labels = risk_2)) %>%  
  mutate(risk = factor(risk, levels = risk_2)) 

# Plot
risk_levels = results_2 %>%
  mutate(risk = factor(p_d_I, labels = risk_2)) %>%  
  mutate(risk = factor(risk, levels = risk_2)) %>% 
  ggplot(aes(x = lambda_FOI*100, y = n)) +
  geom_line(aes(color = factor(risk))) +
  facet_wrap(~ risk, ncol = 1) +
  my_theme +
  labs(
    x = "Attack rate (%)",
    y = "Sample size"
  ) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 2, 2, 2),  # Reduce top margin
    panel.spacing = unit(0.5, "lines")  # Reduce spacing between facets
  ) +
geom_vline(xintercept = 5, linetype = 2) +
geom_hline(
  data = hline_data,
  aes(yintercept = n, group = risk),
  linetype = 2,
  color = "black"
) +
  scale_y_log10() 
# combine and save

out1 = VE_high + risk_levels +
  plot_layout(widths = c(2,1)) +
  plot_annotation(tag_levels = "a")
  
  
# Save
ggsave(out1, 
       filename = "output/sample_size_plot.jpg",
       units = "cm",
       height = 10, width = 14)


# specific events 

# CHIK young 30%
calculate_sample_size_safety(
  p_d_I = all_risk$ifr[1],
  lambda_FOI = c(0.3),
  VE = 0.95,
  max_n = 400000
)

# CHIK old 30%
calculate_sample_size_safety(
  p_d_I = all_risk$ifr[2],
  lambda_FOI = c(0.3),
  VE = 0.95,
  max_n = 400000
)



# SARS
calculate_sample_size_safety(
  p_d_I = 0.01,
  lambda_FOI = c(0.05),
  VE = 0.95,
  max_n = 400000
)


# ebola 1%
calculate_sample_size_safety(
  p_d_I = 0.5,
  lambda_FOI = c(0.05),
  VE = 0.95,
  max_n = 400000
)
