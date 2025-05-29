### Calculate power to detect adverse events in trial 

library(tidyverse)
library(patchwork)
library(stats)

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

all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65-69", "70-79", "80+"))) %>% 
  arrange(age_group)


# Function to calculate sample size for vaccine safety test
calculate_sample_size_safety = function(p_d_I,
                                        lambda_FOI,
                                        VE = 0.95,
                                        criterion = 1.0,
                                        # Additional safety margin
                                        alpha = 0.05,
                                        power = 0.80,
                                        p_adverse_true = NULL) {
  # p0: Maximum acceptable adverse event rate
  p0 = VE * lambda_FOI * p_d_I * criterion
  
  # p1: True adverse event rate we want to detect
  if (is.null(p_adverse_true)) {
    p1 = p0 * 0.5  # Detect if true rate is half the threshold
  } else {
    p1 = p_adverse_true
  }

  # one sided 
  z_alpha = qnorm(1 - alpha)
  z_beta = qnorm(power)
  
  # Sample size formula
  n = (z_alpha * sqrt(p0 * (1 - p0)) + z_beta * sqrt(p1 * (1 - p1)))^2 / 
    (p0 - p1)^2
  
  return(ceiling(n))
}

# Example with realistic parameters
p_d_I = 0.02      # 2% severe disease risk if infected
lambda_FOI = 0.15 # 15% annual infection probability
VE = 0.95         # 95% vaccine efficacy
criterion = 0.5   # Require 2:1 benefit-to-risk ratio 

# Calculate sample size
n = calculate_sample_size_safety(p_d_I, lambda_FOI, VE, criterion)
n

# Create a grid of values
p_d_I_values = seq(0.01,0.2, 0.001)  # Disease severity
lambda_FOI_values = seq(0.05, 0.6, 0.01)  # Infection probability
criterion_values = c(1.0)  # Safety margins

# Calculate for different FOI and criterion combinations
results = expand.grid(
  p_d_I = p_d_I_values, 
  lambda_FOI = lambda_FOI_values,
  criterion = criterion_values
)

results$n = mapply(calculate_sample_size_safety,
                    p_d_I = results$p_d_I,
                    lambda_FOI = results$lambda_FOI,
                    criterion = results$criterion)

results %>%  
  ggplot(aes(x = p_d_I*100000, y = lambda_FOI*100, fill = n)) +  # Add 'fill = n' inside aes()
  geom_tile() +  # Remove 'fill = n' from here
  scale_fill_gradient2(
    low = "#C51B8A", 
    mid = "white",
    high = "#FD8D3C",
    midpoint = 4000) +
  my_theme + 
  labs(x = "Risk from infection (per 100,000)", 
       y = "Attack rate (%)", 
       fill = "Sample Size") 
