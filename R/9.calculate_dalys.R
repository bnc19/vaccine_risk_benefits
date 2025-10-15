#  Script to calculate DALYs averted by IXCHIQ vaccination. 
# File outputs Figure 3C/F.

library(tidyverse)
library(cowplot)

set.seed(12)

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
    strip.text = element_blank()
  )

# Load and prepare data ----
all_risk = read.csv("data/all_risk.csv")[,-1] %>%
  mutate(age_group = factor(age_group, levels = c("18-64", "65+")))

# Constants ----
n_samples = 100000
prob_symp_given_inf = 0.5
prob_chronic_given_sev = 0.5
acute_duration = 7/365
chronic_duration = 1
life_expect = 84

# Disability weights
severe_dw = 0.133
mild_dw = 0.006  
chronic_dw = 0.233

# Vaccination parameters
ve_symp = 0.95 # 0.5
ve_inf = 0 
coverage = 1 # per 10000 vaccinated 
prob_symp_given_inf_vac = (1 - ve_symp) / (1 - ve_inf) * prob_symp_given_inf

# Population 
populations = c(young = 10000, old = 10000) # per 10000 vaccinated 
ages = c(young = 41, old = 75.6)

# Attack rates
endemic_lambda = 0.024
attack_rates = c(
  traveller = 1 - exp(-(-log(0.7)/90) * 10),
  low = 0.05,
  high = 0.3,
  endemic = 1 - (1 - endemic_lambda)^10
)

# Generate probability samples ----
generate_probability_samples = function(risk_data, n_samples) {
  
  # Vectorized sampling for both age groups
  samples = risk_data %>%
    rowwise() %>%
    mutate(
      # Death probabilities
      prob_death_infection = list(rbinom(n_samples, round(incidence), ifr) / round(incidence)),
      prob_death_vaccine = list(rbinom(n_samples, round(doses), p_deaths) / round(doses)),
      
      # Case probabilities  
      prob_case_infection = list(rbinom(n_samples, round(incidence), p_case) / round(incidence)),
      prob_case_vaccine = list(rbinom(n_samples, round(doses), p_sae) / round(doses))
    ) %>%
    ungroup()
  
  return(samples)
}

probability_samples = generate_probability_samples(all_risk, n_samples)

# Extract samples by age group
young_samples = probability_samples[probability_samples$age_group=="18-64", ]
old_samples = probability_samples[probability_samples$age_group=="65+", ]

# Calculate derived probabilities ----
calculate_disease_probabilities = function(prob_case, prob_death) {
  # Ensure probabilities are valid [0,1]
  
  if(any(prob_case > prob_symp_given_inf)) stop("invalid prob")
  
  prob_severe_given_symptomatic = prob_case / prob_symp_given_inf
  prob_mild_given_symptomatic = (1 - prob_severe_given_symptomatic)
  case_fatality_ratio = prob_death / prob_case
  
  return(list(
    prob_severe = prob_severe_given_symptomatic,
    prob_mild = prob_mild_given_symptomatic, 
    cfr = case_fatality_ratio
  ))
}

young_probs = calculate_disease_probabilities(
  young_samples$prob_case_infection[[1]], 
  young_samples$prob_death_infection[[1]]
)

young_probs$prob_sae_vaccine = young_samples$prob_case_vaccine[[1]]
young_probs$prob_death_vaccine = young_samples$prob_death_vaccine[[1]]

old_probs = calculate_disease_probabilities(
  old_samples$prob_case_infection[[1]], 
  old_samples$prob_death_infection[[1]]
)

old_probs$prob_sae_vaccine = old_samples$prob_case_vaccine[[1]]
old_probs$prob_death_vaccine = old_samples$prob_death_vaccine[[1]]


# Vectorized outcome calculation ----
calculate_outcomes_unvac = function(attack_rate, populations, young_probs, old_probs) {
  
  # Calculate infections
  infections = populations * attack_rate
  cases = infections * prob_symp_given_inf
  
  # Vectorized calculations across all samples
  results = tibble(
    sample = 1:n_samples,
    
    # Young population outcomes
    mild_young = cases["young"] * young_probs$prob_mild,
    severe_young = cases["young"] * young_probs$prob_severe,
    chronic_young = severe_young * prob_chronic_given_sev,
    deaths_young = severe_young * young_probs$cfr,
    
    # Old population outcomes  
    mild_old = cases["old"] * old_probs$prob_mild,
    severe_old = cases["old"] * old_probs$prob_severe,
    chronic_old = severe_old * prob_chronic_given_sev,
    deaths_old = severe_old * old_probs$cfr,
    
  )
  
  return(results)
}


calculate_outcomes_vac = function(attack_rate, populations, young_probs, old_probs) {
  
  # Calculate infections
  infections_unvac = populations * attack_rate * (1-coverage)
  cases_unvac = infections_unvac * prob_symp_given_inf
  infections_vac = populations * attack_rate * (1-ve_inf) * coverage
  cases_vac = infections_vac * prob_symp_given_inf_vac
  
  # Vectorized calculations across all samples
  results = tibble(
    sample = 1:n_samples,
    
    # Young population outcomes
    vsae_young = populations["young"] * coverage * young_probs$prob_sae_vaccine,
    vdeaths_young = populations["young"] * coverage * young_probs$prob_death_vaccine,
    
    mild_young = (cases_vac["young"] + cases_unvac["young"]) * young_probs$prob_mild,
    severe_young = (cases_vac["young"] + cases_unvac["young"]) * young_probs$prob_severe,
    chronic_young = severe_young * prob_chronic_given_sev,
    deaths_young = severe_young * young_probs$cfr,
    
    # Old population outcomes  
    vsae_old = populations["old"] * coverage *  old_probs$prob_sae_vaccine,
    vdeaths_old = populations["old"] * coverage * old_probs$prob_death_vaccine,
    
    mild_old = (cases_vac["old"] + cases_unvac["old"]) * old_probs$prob_mild,
    severe_old = (cases_vac["old"] + cases_unvac["old"]) * old_probs$prob_severe,
    chronic_old = severe_old * prob_chronic_given_sev,
    deaths_old = severe_old * old_probs$cfr
  )
  
  return(results)
}

# Calculate results for all attack rates ----
nv_scenario_results = map(attack_rates, ~calculate_outcomes_unvac(.x, populations, young_probs, old_probs))
names(nv_scenario_results) = names(attack_rates)

v_scenario_results = map(attack_rates, ~calculate_outcomes_vac(.x, populations, young_probs, old_probs))
names(v_scenario_results) = names(attack_rates)


# Calculate DALYs ----
calculate_dalys = function(results_df, populations, coverage, v = F) {
  
  
  number_vacc = populations * coverage
  
  # DALYs due to infection 
  results_df = results_df %>%
    mutate(
      # Years of Life Lost (YLL)
      yll_young = deaths_young * (life_expect - ages["young"]),
      yll_old = deaths_old * (life_expect - ages["old"]),
      
      # Years Lived with Disability (YLD)  
      yld_mild_young = mild_young * mild_dw * acute_duration,
      yld_severe_young = severe_young * severe_dw * acute_duration,
      yld_chronic_young = chronic_young * chronic_dw * chronic_duration,
      
      yld_mild_old = mild_old * mild_dw * acute_duration,
      yld_severe_old = severe_old * severe_dw * acute_duration,
      yld_chronic_old = chronic_old * chronic_dw * chronic_duration,
      
      # Total DALYs
      dalys_young = (yll_young + yld_mild_young + yld_severe_young + yld_chronic_young) , 
      dalys_old = (yll_old +  yld_mild_old + yld_severe_old + yld_chronic_old) 
      
    )
  
  if(v==T){
    results_df = results_df %>%
      mutate(
        vyll_young = vdeaths_young * (life_expect - ages["young"]),
        vyll_old = vdeaths_old * (life_expect - ages["old"]), 
        vyld_severe_young = vsae_young * severe_dw * acute_duration,
        vyld_severe_old = vsae_old * severe_dw * acute_duration,
        vdalys_young = (vyll_young + vyld_severe_young) , 
        vdalys_old = (vyll_old + vyld_severe_old) ,
        totaldalys_young = dalys_young + vdalys_young,
        totaldalys_old = dalys_old + vdalys_old
      )
    
  }
  
  return(results_df)
  
}

# Apply DALY calculation to all scenarios
nv_results = map(nv_scenario_results, calculate_dalys, populations = populations, coverage = coverage)
v_results = map(v_scenario_results, calculate_dalys,populations = populations, coverage = coverage, v=T)


calculate_dalys_averted = function(nv_results, v_results){
  
  nv_results$dalys_averted_young = nv_results$dalys_young - v_results$totaldalys_young
  nv_results$dalys_averted_old =  nv_results$dalys_old - v_results$totaldalys_old
  
  return(nv_results)
}

all_results = map2(nv_results, v_results, calculate_dalys_averted)

summarise_results = function(all_results, exclude_cols = NULL) {
  map_dfr(all_results, ~{
    # Get all numeric columns
    numeric_cols <- .x %>% 
      select(where(is.numeric)) %>% 
      names()
    
    # Remove excluded columns if specified
    if (!is.null(exclude_cols)) {
      numeric_cols <- setdiff(numeric_cols, exclude_cols)
    }
    
    # Create summary for each numeric column
    summary_stats <- map_dfc(numeric_cols, function(col) {
      tibble(
        !!paste0("mean_", col) := mean(.x[[col]], na.rm = TRUE),
        !!paste0("low_", col) := quantile(.x[[col]], 0.025, na.rm = TRUE),
        !!paste0("high_", col) := quantile(.x[[col]], 0.975, na.rm = TRUE)
      )
    })
    
    summary_stats
  }, .id = "scenario")
}

final_results = summarise_results(all_results, exclude_cols = "sample")
vaccine_results = summarise_results(v_results, exclude_cols = "sample")


# Plot dalys without vaccination 

p_left = final_results %>% 
  select(scenario, mean_dalys_young:high_dalys_old) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_replace(name, "_dalys_", ".")) %>% 
  separate(name, into = c("name", "age")) %>%  
  mutate(value =  - value) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(scenario = factor(scenario, levels = c("low","high","endemic","traveller"),
                           labels = c("Small\noutbreak",
                                      "Large\noutbreak",
                                      "Endemic",
                                      "Traveller")),
         age = factor(age, levels = c("young", "old"), labels = c("18-64", "65+"))) %>% 
  ggplot(aes(y = mean, x = age)) +
  geom_bar(stat = "identity", fill = "#A0BFD8") +
  geom_errorbar(aes(ymin = low, ymax = high), width=0.25) +
  coord_flip(clip = "off") +
  facet_wrap(~ scenario, ncol = 1) +
  my_theme  +
  ylim(c(-100,0))

# Plot Dalys with vaccination 

vaccine_dalys = vaccine_results %>% 
  select(scenario, contains("dalys")) %>%  
  pivot_longer(-scenario) %>% 
  separate(name, into = c("name", "outcome", "age")) %>%  
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(scenario = factor(scenario, levels = c("low","high","endemic","traveller"),
                           labels = c("Small\noutbreak",
                                      "Large\noutbreak",
                                      "Endemic",
                                      "Traveller")),
         age = factor(age, levels = c("young", "old"), labels = c("18-64", "65+"))) 

# Filter for mean deaths attributal to vaccine or infection 
vdaly_mean_data = vaccine_dalys %>% 
  filter(outcome %in% c("dalys", "vdalys")) %>% 
  select(-low, -high)

# Filter for error bars (from total deaths only (vaccine + infection))
vdaly_error_data = vaccine_dalys %>% 
  filter(outcome == "totaldalys") %>% 
  select(-mean, -outcome) # - outcome so doesn't join by outcome 

# Merge error bars into mean data
vdaly_plot_data = vdaly_mean_data %>% 
  left_join(vdaly_error_data)


p_right = vdaly_plot_data %>% 
  ggplot(aes(y = mean, x = age, fill = outcome)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low, ymax = high), width=0.25) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("vdalys" = "#A3B18A",
                               "dalys" = "#4B4B4B")) +
  my_theme +
  facet_wrap(~ scenario, ncol = 1) +
  labs(y = "DALYs (per 10,000)") +
  ylim(c(0,100))

# Plot DALYS averted 

dalys_averted = final_results %>%  
  select(scenario, mean_dalys_averted_young:high_dalys_averted_old) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_replace(name, "_dalys_averted_", ".")) %>% 
  separate(name, into = c("name", "age")) %>%  
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(scenario = factor(scenario, levels = c("low","high","endemic","traveller"),
                           labels = c("Small\noutbreak",
                                      "Large\noutbreak",
                                      "Endemic",
                                      "Traveller")),
         age = factor(age, levels = c("young", "old"), labels = c("18-64", "65+"))) %>% 
  mutate(color = ifelse(mean <0, "0", "1")) %>% 
  ggplot(aes(x = mean, y = age)) +
  geom_bar(stat = "identity", aes(fill = color)) +
  geom_errorbar(aes(xmin = low, xmax = high),  width =0.25) +
  facet_wrap(~ scenario ,  ncol = 1) +
  my_theme +
  scale_fill_manual(values= c("#1B3B6F")) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlim(c(-100,100))



ggsave(dalys_averted, 
       filename = "output/difference_in_dalys.pdf",
       units = "cm",
       height = 8,
       width = 10,
)


combined_plot = plot_grid(p_left, p_right)

ggsave(combined_plot,
       filename = "output/dalys.pdf",
       units = "cm",
       height = 8,
       width = 10)

