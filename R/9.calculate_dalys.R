# Script to calculate DALYs averted by IXCHIQ vaccination. 
# File outputs Figure 3C/F.

library(tidyverse)

set.seed(12)

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
    legend.title = element_text(size = 7)
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
coverage = 0.25
prob_symp_given_inf_vac = (1 - ve_symp) / (1 - ve_inf) * prob_symp_given_inf

# Population 
populations = c(young = 489834, old = 110157) 
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
    
    # So names match vaccine scenario
    deaths_young_total = deaths_young,
    deaths_old_total = deaths_old,
    severe_young_total = severe_young,
    severe_old_total = severe_old
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
    
    severe_young_total = severe_young + vsae_young,
    deaths_young_total = deaths_young + vdeaths_young,
    
    # Old population outcomes  
    vsae_old = populations["old"] * coverage *  old_probs$prob_sae_vaccine,
    vdeaths_old = populations["old"] * coverage * old_probs$prob_death_vaccine,
    
    mild_old = (cases_vac["old"] + cases_unvac["old"]) * old_probs$prob_mild,
    severe_old = (cases_vac["old"] + cases_unvac["old"]) * old_probs$prob_severe,
    chronic_old = severe_old * prob_chronic_given_sev,
    deaths_old = severe_old * old_probs$cfr,
    
    severe_old_total = severe_old + vsae_old, # severe disease and sae are equivalent 
    deaths_old_total = deaths_old + vdeaths_old # deaths are equivalent 
  )
  
  return(results)
}

# Calculate results for all attack rates ----
nv_scenario_results = map(attack_rates, ~calculate_outcomes_unvac(.x, populations, young_probs, old_probs))
names(nv_scenario_results) = names(attack_rates)

v_scenario_results = map(attack_rates, ~calculate_outcomes_vac(.x, populations, young_probs, old_probs))
names(v_scenario_results) = names(attack_rates)


# Calculate DALYs ----
calculate_dalys = function(results_df) {
  results_df %>%
    mutate(
      # Years of Life Lost (YLL)
      yll_young = deaths_young_total * (life_expect - ages["young"]),
      yll_old = deaths_old_total * (life_expect - ages["old"]),
      
      # Years Lived with Disability (YLD)  
      yld_mild_young = mild_young * mild_dw * acute_duration,
      yld_severe_young = severe_young_total * severe_dw * acute_duration,
      yld_chronic_young = chronic_young * chronic_dw * chronic_duration,
      
      yld_mild_old = mild_old * mild_dw * acute_duration,
      yld_severe_old = severe_old_total * severe_dw * acute_duration,
      yld_chronic_old = chronic_old * chronic_dw * chronic_duration,
      
      # Total DALYs
      dalys_young = yll_young + yld_mild_young + yld_severe_young + yld_chronic_young,
      dalys_old = yll_old +  yld_mild_old + yld_severe_old + yld_chronic_old
      
    )
}

# Apply DALY calculation to all scenarios
nv_results = map(nv_scenario_results, calculate_dalys)
v_results = map(v_scenario_results, calculate_dalys)

calculate_dalys_averted = function(nv_results, v_results, populations, coverage){
  
  nv_results$dalys_averted_young = nv_results$dalys_young - v_results$dalys_young
  nv_results$dalys_averted_old = nv_results$dalys_old - v_results$dalys_old
  
  number_vacc = populations * coverage
  nv_results$dalys_averted_young_per_10000 =  nv_results$dalys_averted_young / number_vacc["young"] * 10000
  nv_results$dalys_averted_old_per_10000 =  nv_results$dalys_averted_old / number_vacc["old"] * 10000
  
  
  return(nv_results)
}

all_results = map2(nv_results, v_results, calculate_dalys_averted, populations = populations, coverage = coverage)


# Summary statistics ----
summarise_results = function(all_results) {
  map_dfr(all_results, ~{
    .x %>%
      summarise(
        mean_young = mean(dalys_averted_young),
        low_young = quantile(dalys_averted_young, 0.025),
        high_young = quantile(dalys_averted_young, 0.975),
        mean_old = mean(dalys_averted_old),
        low_old = quantile(dalys_averted_old, 0.025),
        high_old = quantile(dalys_averted_old, 0.975),
        mean_young2 = mean(dalys_averted_young_per_10000),
        low_young2 = quantile(dalys_averted_young_per_10000, 0.025),
        high_young2 = quantile(dalys_averted_young_per_10000, 0.975),
        mean_old2 = mean(dalys_averted_old_per_10000),
        low_old2 = quantile(dalys_averted_old_per_10000, 0.025),
        high_old2 = quantile(dalys_averted_old_per_10000, 0.975),
        .groups = "drop"
      )
  }, .id = "scenario")
}

final_results = summarise_results(all_results)
print(final_results)


# Plot DALYS averted 

final_plot = final_results %>%  
  select(-c(mean_young:high_old)) %>% 
  pivot_longer(-scenario) %>% 
  separate(name, into = c("name", "age")) %>%  
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(scenario = factor(scenario, levels = c("low","high","endemic","traveller"),
                           labels = c("Small\noutbreak",
                                      "Large\noutbreak",
                                      "Endemic",
                                      "Traveller")),
         age = factor(age, levels = c("young2", "old2"), labels = c("18-64", "65+"))) %>% 
  ggplot(aes(x = scenario, y = mean)) +
  geom_col(aes(fill = scenario)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.1)+
  my_theme +
  labs(y="DALYs averted per 10,000 vaccinated", x = "") +
  scale_fill_manual(values = c("#E6A23C", "#E06685", "#9689C3", "#5BB1C2")) +
  theme( legend.position = "none") +
  facet_wrap(~age) +
  geom_hline(yintercept = 0, linetype = 2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.4) +
  scale_y_continuous(breaks = c(-20, 0, 20, 40 , 60, 80))


if(ve_symp != 0.5){
ggsave(final_plot, 
       filename = "output/plot_dalys.jpg",
       units = "cm",
       height = 5,
       width = 12)
} else{
  ggsave(final_plot, 
         filename = "output/plot_dalys_ve50.jpg",
         units = "cm",
         height = 5,
         width = 12) 
}
final_plot
