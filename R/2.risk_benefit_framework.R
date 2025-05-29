library(tidyverse)
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

# Calculate death risk
# --- import data
all_risk = read.csv("data/all_risk.csv")[,-1]
all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) %>% 
  arrange(age_group)

# define baseline parameters
VE = c(0.5, 0.95)

endemic_lambda = 0.024 # Gabriel's paper
endemic_AR = (1-(1-endemic_lambda)^10)
low_AR = 0.05
high_AR = 0.3
days_hol = 10
# traveller in high incidence outbreak
# 0.3 = 1- exp(-90 * r)
r = - log(0.7)/90
traveller_AR = 1  - exp(-r * days_hol)
AR = c(traveller_AR, low_AR, high_AR,  endemic_AR)

population = 10000

################################################################################
################# Infection / Vaccine risk sensitivity analysis ################
################################################################################
# see what chikv range is 
all_risk$p_case * population
all_risk$p_deaths  * population
all_risk$p_sae  * population
all_risk$ifr  * population

# Define scaling factors for sensitivity analysis
inf_cases = seq(0,2000,by=50)
vac_cases = seq(0,20,by= 1)

# Create an empty list to store all results
sensitivity_results_cases = sensitivity_results_deaths = list()
scenario_counter = 1

for (ve_idx in 1:length(VE)) {
  
  current_VE = VE[ve_idx]
  
  for (ar_idx in 1:length(AR)) {

    current_AR = AR[ar_idx]
    
  
  # Loop through sensitivity parameters
  for (ic in inf_cases) {
      for (vc in vac_cases) {
          # Create a copy of the original dataframe for this iteration
          temp_risk = all_risk
          
          
          # Calculate risk metrics using the current AR value and scaled parameters
          temp_risk$risk_dis_nv = ic * current_AR  
          temp_risk$risk_dis_v = ic * current_AR * (1-current_VE) + vc
          temp_risk$cases_averted = temp_risk$risk_dis_nv - temp_risk$risk_dis_v
          
          # Add parameters and scaling information
          temp_risk$daily_ar = current_AR
          temp_risk$ve = current_VE
          temp_risk$ar_scenario = ar_idx
          temp_risk$inf_cases = ic
          temp_risk$vac_cases = vc
          temp_risk$scenario_id = scenario_counter
          
          # Store the results for this iteration
          sensitivity_results_cases[[scenario_counter]] = temp_risk
          scenario_counter = scenario_counter + 1
        }
      }
    }
}

# Combine all results
c_sensitivity_results = do.call(rbind, sensitivity_results_cases)

c_sensitivity_results = c_sensitivity_results %>%  
  mutate(ar_scenario = factor(ar_scenario, levels = c(2,3,4,1),
                              labels = c("Small outbreak",
                                         "Large outbreak",
                                         "Endemic",
                                         "Traveller"))) %>% 
  mutate(age_group = fct_rev(age_group),
         ve = factor(ve, levels = c(0.5,0.95), labels = c("50% vaccine efficacy", "95% vaccine efficacy")))



# add chik points 
points_df = bind_rows(
  c_sensitivity_results %>%
    mutate(x_val = p_sae * population,
           y_val = p_case * population,
           source = "SAE"),
  
  c_sensitivity_results %>%
    mutate(x_val = p_deaths * population,
           y_val = ifr * population,
           source = "Deaths")
) %>%
  mutate(point_type = interaction(age_group, source, sep = " "))

c_SA = c_sensitivity_results %>%
  ggplot(aes(x = vac_cases, y = inf_cases)) +
  geom_tile(aes(fill = cases_averted)) +
  facet_grid(ar_scenario ~ ve) +
  scale_fill_gradient2(
    low = "#C51B8A",
    mid = "white",
    high = "#FD8D3C",
    midpoint = 0,
    limits = c(min(c_sensitivity_results$cases_averted), max(c_sensitivity_results$cases_averted)),
    trans = "pseudo_log",
    breaks = c(-50, -5, 0, 5, 50, 500)
  ) +
  geom_contour(aes(z = cases_averted),
               breaks = 0,
               linetype = 2,
               color = "black",
               linewidth = 0.2) +
  geom_point(data = points_df, 
             aes(x = x_val, y = y_val, shape = point_type),
             color = "black", size = 0.8, stroke = 0.2) +
  labs(x = "Risk of vaccine SAE (per 10,000)",
       y = "Risk of severe outcomes following infection (per 10,000)",
       fill = "Outcomes averted\n(per 10,000)",
       shape = "Chikungunya\nadverse events") +
  my_theme +
  scale_shape_manual(values = c(0,15,2,17)) +
  guides(shape = guide_legend(override.aes = list(size = 2)))

c_SA

ggsave(c_SA,
       filename = paste0("output/cases_averted_plot.jpeg"),
       units = "cm",
       height = 10, width = 12)
