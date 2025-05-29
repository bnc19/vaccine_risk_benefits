#### Calculate probability benefits greater than risk for chik trial 

library(tidyverse)
library(scales)  # for comma()
library(patchwork)

cols = c(  "#FD8D3C", 
           "#F768A1",
           "#C51B8A",
           "#7A0177" )


# --- import data

all_risk = read.csv("data/all_risk.csv")[,-1]

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

all_risk = all_risk %>%  
  mutate(age_group = factor(age_group, levels = c("18-64", "65+"))) %>% 
  mutate(age_group = fct_rev(age_group))

# define baseline parameters

endemic_lambda = 0.024 # Gabriel's paper
endemic_AR = (1-(1-endemic_lambda)^10)
low_AR = 0.05
high_AR = 0.3
days_hol = 10
# traveller in Large outbreak
# 0.3 = 1- exp(-09 * r)
r = - log(0.7)/90
traveller_AR = 1  - exp(-r * days_hol)
AR = c(traveller_AR, low_AR, high_AR,  endemic_AR)


# define observed risks 
all_risk$trial_pop =  c(2736, 346)# from table 1 https://www.sciencedirect.com/science/article/pii/S0140673623006414?via%3Dihub#sec1
 

### Function to calculate probability that pv > p1 given observed outcomes and n sample size 
calculate_prob_bayesian = function(x, n, p1) {
  prob = pbeta(p1, x + 0.5, n - x + 0.5)  
  return(prob)
}


# Function to calculate probabilities for a specific outcome
calculate_outcome_probs = function(n_vac_outcome, p_outcome_inf, outcome_name, n, AR, VE = 0.95) {
  
  p_outcome_given_risk = data.frame()
  
  
  for (k in seq_along(AR)) {
    p_row = numeric(nrow(all_risk))
    for (i in 1:nrow(all_risk)) {
      p_row[i] = calculate_prob_bayesian(
        x = n_vac_outcome[i],  # Use the outcome column (v_deaths or v_sae)
        n = n[i],
        p1 = p_outcome_inf[i] * VE * AR[k]  # threshold risk 
        
      )
    }
    p_outcome_given_risk = rbind(p_outcome_given_risk, 
                                 data.frame(AR = AR[k], 
                                            outcome = outcome_name,
                                            scenario = k, 
                                            t(p_row)))
  }
  
  # Give column names to the results
  colnames(p_outcome_given_risk) = c("AR", "outcome", "scenario", paste(all_risk$age_group))
  return(p_outcome_given_risk)
}


################################ Generalized plot ##############################

p1 = expand_grid(ar = AR, 
                 VE = c(0.5, 0.95), 
                 pop = c(1000, 10000), 
                 p_event = c(0, 0.001, 0.01), 
                 inf_risk = seq(0, 0.3, 0.01)) %>% 
  mutate(p1 = inf_risk * VE * ar) %>% 
  mutate(p_benefit = calculate_prob_bayesian(p_event*pop, pop, p1)) %>% 
  mutate(pop = paste0(comma(pop))) %>%  # add comma here
  mutate(ar = factor(ar, levels = c(AR[2],AR[3],AR[4],AR[1]),
                              labels = c("Small outbreak",
                                         "Large outbreak",
                                         "Endemic",
                                         "Traveller"))) %>% 
  ggplot(aes(y = p_benefit, x = inf_risk)) +
  geom_line(aes(color = factor(p_event), linetype = pop), size =  0.7) +
  facet_grid(paste0(VE*100, "% vaccine efficacy") ~ ar) +
  labs(x = "Probability of severe outcomes following infection", 
       y = "Probability vaccine benefit",
       color = "Probability of\nvaccine SAE",
       linetype = "Cohort size") +
  my_theme +
  theme(legend.position = "top") + 
  scale_color_manual(values = c("#74C476", "#6A51A3",  "#67A9CF"))

p1


ggsave(
  p1,
  filename = "output/prob_vac_benefit_general.jpeg",
  units = "cm",
  height = 9,
  width = 14
)


########################### chikungunya specific ###############################

VE = 0.95 # vary here for 0.5 VE

# Calculate for deaths 
p_deaths = calculate_outcome_probs(n_vac_outcome = all_risk$v_deaths,
                                   p_outcome_inf = all_risk$ifr,
                                   outcome_name = "deaths_All data",
                                   n = all_risk$doses,
                                   AR, 
                                   VE = VE)

p_deaths_trial = calculate_outcome_probs(c(0,0),
                                         p_outcome_inf = all_risk$ifr,
                                         outcome_name = "deaths_Trial data",
                                         n = all_risk$trial_pop,
                                         AR, 
                                         VE = VE
)


# Calculate for SAEs 
p_sae = calculate_outcome_probs(n_vac_outcome = all_risk$SAE,
                                p_outcome_inf = all_risk$p_case,
                                outcome_name = "SAE_All data",
                                n = all_risk$doses,
                                AR, 
                                VE = VE)

p_sae_trial = calculate_outcome_probs(n_vac_outcome = c(1,1),
                                p_outcome_inf = all_risk$p_case,
                                outcome_name = "SAE_Trial data",
                                n = all_risk$trial_pop, 
                                AR, 
                                VE = VE)


# Combine the results
p_combined = rbind(p_deaths, p_sae, p_deaths_trial, p_sae_trial)

data_plot = p_combined %>% 
  pivot_longer(-c(AR, scenario, outcome), names_to = "age_group") %>% 
  mutate(scenario = factor(scenario, levels = c(2,3,4,1),
                           labels = c("Small outbreak",
                                      "Large outbreak",
                                      "Endemic",
                                      "Traveller"))) %>%  
  separate(outcome, into= c("outcome", "observed"), sep = "_") %>%  
  mutate(outcome = factor(outcome, labels = c("Deaths", "SAE/clinically-attended chikungunya"))) %>% 
  mutate(scenario = fct_rev(scenario),
         outcome = fct_rev(outcome)) 


p2 = data_plot %>% 
  filter(outcome != "Deaths") %>% 
ggplot(aes(x = value, y = age_group)) +
  geom_point(aes(color = scenario, shape = observed), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", x = "Probability vaccine benefit\noutweighs risk of SAE") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  theme(legend.position = "none")


p3 = data_plot %>% 
  filter(outcome == "Deaths") %>% 
  ggplot(aes(x = value, y = age_group)) +
  geom_point(aes(color = scenario, shape = observed), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(x = "Probability vaccine benefit\noutweighs risk of death") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank()
  )

out = p2+p3+ 
  plot_annotation(tag_levels = "a")


ggsave(
  out,
  filename = paste0("output/prob_vac_benefit_VE_", VE, ".jpeg"),
  units = "cm",
  height = 6,
  width = 14
)



# sensitivity analysis with 50% lower p_case -----------------------------------

# Calculate for SAEs 
p_sae_2 = calculate_outcome_probs(n_vac_outcome = all_risk$SAE,
                                p_outcome_inf = all_risk$p_case/2,
                                outcome_name = "SAE_All data",
                                n = all_risk$doses,
                                AR, 
                                VE = 0.95)

p_sae_trial_2 = calculate_outcome_probs(n_vac_outcome = c(1,1),
                                      p_outcome_inf = all_risk$p_case/2,
                                      outcome_name = "SAE_Trial data",
                                      n = all_risk$trial_pop, 
                                      AR, 
                                      VE = 0.95)



# Combine the results
p_sa = rbind(p_sae_2,  p_sae_trial_2)

data_plot_sa = p_sa %>% 
  pivot_longer(-c(AR, scenario, outcome), names_to = "age_group") %>% 
  mutate(scenario = factor(scenario, levels = c(2,3,4,1),
                           labels = c("Small outbreak",
                                      "Large outbreak",
                                      "Endemic",
                                      "Traveller"))) %>%  
  separate(outcome, into= c("outcome", "observed"), sep = "_") %>%  
  mutate(scenario = fct_rev(scenario)) 


p4 = data_plot_sa %>% 
  filter(outcome != "Deaths") %>% 
  ggplot(aes(x = value, y = age_group)) +
  geom_point(aes(color = scenario, shape = observed), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", x = "Probability vaccine benefit\noutweighs risk of SAE") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  theme(legend.title = element_blank())



ggsave(
  p4,
  filename = "output/prob_vac_benefit_half_p_case.jpeg",
  units = "cm",
  height = 6,
  width = 14
)




# sensitivity analysis with 1 death --------------------------------------------

p_deaths_sa = calculate_outcome_probs(n_vac_outcome = c(0,1),
                                   p_outcome_inf = all_risk$ifr,
                                   outcome_name = "1 death 65+ years",
                                   n = all_risk$doses,
                                   AR, 
                                   VE = 0.95)


p_deaths_sa2 = calculate_outcome_probs(n_vac_outcome = c(0,0),
                                      p_outcome_inf = all_risk$ifr,
                                      outcome_name = "0 deaths 65+ years",
                                      n = all_risk$doses,
                                      AR, 
                                      VE = 0.95)




data_plot_sa_deaths = p_deaths_sa %>% 
  bind_rows(p_deaths_sa2) %>%  
  pivot_longer(-c(AR, scenario, outcome), names_to = "age_group") %>% 
  mutate(scenario = factor(scenario, levels = c(2,3,4,1),
                           labels = c("Small outbreak",
                                      "Large outbreak",
                                      "Endemic",
                                      "Traveller"))) %>%  
  mutate(scenario = fct_rev(scenario)) 


p5 = data_plot_sa_deaths %>% 
  ggplot(aes(x = value, y = age_group)) +
  geom_point(aes(color = scenario, shape = outcome), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", 
       x = "Probability vaccine benefit\noutweighs risk of death") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  theme(legend.title = element_blank())


ggsave(
  p5,
  filename = "output/prob_vac_benefit_1_death.jpeg",
  units = "cm",
  height = 6,
  width = 14
)


# what if they had recruited more >65 individuals

p_sae_trial_size = calculate_outcome_probs(n_vac_outcome = c(4,1),
                                      p_outcome_inf = all_risk$p_case,
                                      outcome_name = "SAE_Trial data",
                                      n =  c(4000, 1000), 
                                      AR, 
                                      VE = 0.95)


