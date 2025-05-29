library(tidyverse)
library(patchwork)

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
# see what current range is 
all_risk$p_case * population
all_risk$p_deaths  * population
all_risk$p_sae  * population
all_risk$ifr  * population

# Define scaling factors for sensitivity analysis
inf_cases = seq(0,3000,by=50)
vac_cases = seq(0,50,by= 1)

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
         ve = factor(ve, levels = c(0.5,0.95), labels = c("95% vaccine efficacy", "50% vaccine efficacy"))) %>% 
  mutate(ve = fct_rev(ve))



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


################################################################################
########################### VE sensitivity analysis ############################
################################################################################

# define baseline parameters
VE_sa = seq(0,1,0.05)

# Create an empty list to store all results
sensitivity_results_VE = list()
scenario_counter = 1

# Loop through each AR scenario
for (ar_idx in 1:length(AR)) {
  current_AR = AR[ar_idx]
  
  # Loop through sensitivity parameters
 for (temp_VE in VE_sa){

       # Create a copy of the original dataframe for this iteration
          temp_risk = all_risk
          
          temp_risk$risk_death_nv = temp_risk$ifr * current_AR * population
          temp_risk$risk_death_v = (temp_risk$ifr * current_AR  * (1-temp_VE) + temp_risk$p_deaths) * population
          temp_risk$risk_dis_nv = temp_risk$p_case * current_AR * population 
          temp_risk$risk_dis_v = (temp_risk$p_case * current_AR  * (1-temp_VE) + temp_risk$p_sae) * population
          

          temp_risk$deaths_averted = temp_risk$risk_death_nv - temp_risk$risk_death_v
          temp_risk$cases_averted = temp_risk$risk_dis_nv - temp_risk$risk_dis_v
          
          # Add parameters and scaling information
          temp_risk$daily_ar = current_AR
          temp_risk$ar_scenario = ar_idx
          temp_risk$VE = temp_VE
          temp_risk$scenario_id = scenario_counter
          
          # Store the results for this iteration
          sensitivity_results_VE[[scenario_counter]] = temp_risk
          scenario_counter = scenario_counter + 1
        
          }
      }


# Combine all results
all_sensitivity_results_VE = do.call(rbind, sensitivity_results_VE)


all_sensitivity_results_VE = all_sensitivity_results_VE %>%  
  mutate(ar_scenario = factor(ar_scenario, levels = c(2,3,4,1),
                              labels = c("Small outbreak",
                                         "Large outbreak",
                                         "Endemic",
                                         "Traveller"))) %>% 
  mutate(age_group = fct_rev(age_group)) %>% 
  arrange(age_group)

SA_VE = all_sensitivity_results_VE %>% 
  pivot_longer(cols = c(deaths_averted, cases_averted)) %>% 
  mutate(name = factor(name, levels = c("cases_averted", "deaths_averted"),
                       labels = c("Clinically-attented\nchikungunya cases\n/ SAEs", "Deaths")))

SA_VE_p = SA_VE %>% 
  ggplot(aes(x = VE*100, y = value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group), size = .8) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 95, linetype = 2, colour = "black") +
  labs(x = "Vaccine efficacy (%)", 
       y = "Cases averted (per 10,000 vaccinated)") +
  my_theme +
  theme(legend.position = c(0.08,0.85),
        legend.title = element_blank()) +
  scale_x_continuous(limits  = c(0,100))+
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#F768A1",
    "#C51B8A",
    "#7A0177"  
  )) +
  facet_grid(name~ar_scenario, scales = "free_y")


ggsave(SA_VE_p, 
       filename = "output/SA_VE.jpeg",
       units = "cm",
       height = 7,
       width = 14)


################################################################################
########################## FOI sensitivity analysis ############################
################################################################################

AR_SA = seq(0.01,0.5,0.01)

sensitivity_results_3 = list()

for (i in 1:length(AR_SA)) {
  
  
  # Create a copy of the original dataframe for this iteration
  temp_risk = all_risk
  
  # Calculate risk metrics using the current AR value
  temp_risk$risk_death_nv = temp_risk$ifr * AR_SA[i] * population
  temp_risk$risk_death_v = (temp_risk$ifr * AR_SA[i]  * (1-VE) + temp_risk$p_deaths) * population
  temp_risk$risk_dis_nv = temp_risk$p_case * AR_SA[i] * population 
  temp_risk$risk_dis_v = (temp_risk$p_case * AR_SA[i]  * (1-VE) + temp_risk$p_sae) * population
  
  temp_risk$deaths_averted = temp_risk$risk_death_nv - temp_risk$risk_death_v
  temp_risk$cases_averted = temp_risk$risk_dis_nv - temp_risk$risk_dis_v
  
  # Add AR value as a column for reference
  temp_risk$daily_ar = AR_SA[i]

    # Store the results for this iteration
  sensitivity_results_3[[i]] = temp_risk
}

# Combine all results
all_sensitivity_results_AR = do.call(rbind, sensitivity_results_3)

all_sensitivity_results_AR = all_sensitivity_results_AR %>%  
  mutate(age_group = fct_rev(age_group)) %>% 
  arrange(age_group)

ar_df = data.frame(
  ar_values = AR*100,
  # Create two rows for each AR value - one for each facet
  name = rep(factor(c("Severe cases", "Deaths")), each = length(AR))
)

AR_SA_to_plot = all_sensitivity_results_AR %>% 
  pivot_longer(cols = c(deaths_averted, cases_averted)) %>% 
  mutate(name = factor(name, levels = c("cases_averted", "deaths_averted"),
                       labels = c("Severe cases", "Deaths"))) 

AR_p1 = AR_SA_to_plot %>% 
  filter(name != "Deaths") %>% 
  ggplot(aes(x = daily_ar*100, y = value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  mapply(function(ar_val) {
    list(
      annotate("segment", 
               x = ar_val*100, xend = ar_val*100,
               y = -Inf, yend = -20,
               arrow = arrow(length = unit(0.1, "cm")))
    )
  }, AR) +
  labs(x = "",
       y = "Clinically-attented chikungunya\ncases / SAEs averted\n(per 10,000 vaccinated)") +
  my_theme +
  theme(legend.position = c(0.09,0.8),
        legend.title = element_blank()) +
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#F768A1",
    "#C51B8A",
    "#7A0177"  
  )) +
  scale_y_continuous(breaks = c(0, 200, 400, 600))

AR_p2 = AR_SA_to_plot %>% 
  filter(name == "Deaths") %>% 
  ggplot(aes(x = daily_ar*100, y = value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = age_group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  mapply(function(ar_val) {
    list(
      annotate("segment", 
               x = ar_val*100, xend = ar_val*100,
               y = -Inf, yend = -2,
               arrow = arrow(length = unit(0.1, "cm")))
    )
  }, AR) +
  labs(x = "Attack rate (%)",
       y = "Deaths averted\n(per 10,000 vaccinated)") +
  my_theme +
  theme(legend.position = "none") +
  scale_color_manual(values = c(
    "#FD8D3C", 
    "#F768A1",
    "#C51B8A",
    "#7A0177"  
  )) +
  scale_y_continuous(breaks = c(-2, 0, 2, 4))

SA_ar = AR_p1/AR_p2

SA_ar

ggsave(AR_p1/AR_p2, 
       filename = "output/AR_SA_plot.jpg",
       units = "cm",
       height = 10,
       width = 14)
  