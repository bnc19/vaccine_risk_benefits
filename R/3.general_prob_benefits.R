
library(Hmisc)
library(truncnorm)
library(tidyverse)
library(cowplot)
library(patchwork)

source("R/functions.R")

dir.create("output", showWarnings = FALSE)

my_theme = theme_bw() +
  theme(
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.spacing.y = unit(0, "cm"),  # Reduced to zero
    legend.spacing.x = unit(0, "cm"),  # Added to control horizontal spacing
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "right")


n_sim = 100000  # posterior samples
pop = 10000 # cases averted per pop vaccinated 

# define baseline parameters
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
VE = c(0.5, 0.95)

# Extract beta parameters
x = c(0, 0.001, 0.01, 0, 0.001, 0.01)
n = c(1000,1000,1000,10000,10000,10000)
d = seq(0,0.3,0.01)

ci = extract_beta_params(x = d*10000, n = rep(10000,length(d)))
cv = extract_beta_params(x = x*n, n = n)

# Initialize storage for results

summary_results <- tibble()

# Nested loops over i, j, and k
for(h in 1:length(d)){
  for (i in 1:length(VE)) {
    for (j in 1:length(n)) {
      for (k in 1:length(AR)) {
    
      
      # Generate random samples
      pci = rbeta(n_sim, shape1 = ci$alpha[h], shape2 = ci$beta[h])
      pcv = rbeta(n_sim, shape1 = cv$alpha[j], shape2 = cv$beta[j])
      
      # Calculate threshold and prob benefit
      c_threshold = pci * VE[i] * AR[k]
      c_benefit = c_threshold > pcv
    
      # Calculate proportion averted 
      c_av = (c_threshold - pcv) * pop
      
      # Calculate cases with and without vaccination 
      nci = pci * AR[k] * pop 
      ncv = nci  * (1-VE[i])  

      VE_value = VE[i]
      AR_value = AR[k]
      n_value = paste0(x[j], "_", n[j])
      d_value = d[h]
      # Store summary results
      summary_results <- bind_rows(summary_results, 
                                   tibble(
                                     n_value = n_value,
                                     d_value = d_value,
                                     VE_value = VE_value,
                                     AR_value = AR_value,
                                     c_benefit = mean(c_benefit),
                                     mean_c_av = mean(c_av),
                                     low_c_av  = quantile(c_av, probs = c(0.025)),
                                     high_c_av = quantile(c_av, probs = c(0.975)),
                                     mean_pc_i = mean(pci* pop),
                                     mean_pc_v = mean(pcv* pop),
                                     low_pc_i  = quantile(pci* pop, probs = c(0.025)),
                                     low_pc_v  = quantile(pcv* pop, probs = c(0.025)),
                                     high_pc_i = quantile(pci* pop, probs = c(0.975)),
                                     high_pc_v = quantile(pcv* pop, probs = c(0.975)),
                                     mean_nc_i = mean(nci),
                                     mean_nc_v = mean(ncv),
                                     low_nc_i  = quantile(nci, probs = c(0.025)),
                                     low_nc_v  = quantile(ncv, probs = c(0.025)),
                                     high_nc_i = quantile(nci, probs = c(0.975)),
                                     high_nc_v = quantile(ncv, probs = c(0.975)),
                                   ))
          }
  }
}
}

# check probability of benefits are consistent

p1 = summary_results %>% 
  select(AR_value, n_value, d_value, VE_value, c_benefit) %>% 
  mutate(AR_value = factor(AR_value, levels = c(AR[2],AR[3],AR[4],AR[1]),
                     labels = c("Small outbreak",
                                "Large outbreak",
                                "Endemic",
                                "Traveller"))) %>% 
  separate(n_value, into = c("p_event", "pop"), sep = "_") %>% 
  ggplot(aes(y = c_benefit, x = d_value)) +
  geom_line(aes(color = factor(p_event), linetype = pop), size =  0.7) +
  facet_grid(paste0(VE_value*100, "% vaccine efficacy") ~ AR_value) +
  labs(x = "Probability of severe outcome following infection", 
       y = "Probability vaccine benefit",
       color = "Probability of\nvaccine SAE",
       linetype = "Cohort size") +
  my_theme +
  theme(legend.position = "top") + 
  scale_color_manual(values = c("#74C476", "#6A51A3",  "#67A9CF"))



ggsave(
  p1,
  filename = "output/prob_vac_benefit_general_new.jpeg",
  units = "cm",
  height = 9,
  width = 14
)
