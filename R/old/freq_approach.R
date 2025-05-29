
library(tidyverse)
all_risk = read.csv("data/all_risk.csv")[,-1]
all_risk$trial_pop =  c(2736, 346)
all_risk$trial_deaths =  c(0, 0)
all_risk$trial_sae =  c(1,1)

cols = c(  "#FD8D3C", 
           "#F768A1",
           "#C51B8A",
           "#7A0177" )


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



VE = 0.95

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




calculate_prob_freq = function(xi,ni,xv, nv){
  
  prob = prop.test(c(xi, xv), c(ni, nv), alternative = "greater", correct = FALSE)
  
  return(prob$p.value)
}


calculate_prob_freq(xi = all_risk$deaths[2] * 0.03885553 * 00.95, ni = all_risk$incidence[2], xv = 0, nv = all_risk$trial_pop[2])


all_results = data.frame(
  d = as.numeric(),
  c = as.numeric(),
  age = as.character(),
  ar = as.integer()
)

for(j in 1:2) {
  for (i in 1:length(AR)) {
    d = calculate_prob_freq(
      xi = all_risk$deaths[j] * VE * AR[i],
      xv = all_risk$v_deaths[j],
      nv = all_risk$doses[j],
      ni = all_risk$incidence[j]
    )
    
    c = calculate_prob_freq(
      xi = all_risk$cases[j] * VE * AR[i],
      xv = all_risk$SAE[j],
      nv = all_risk$doses[j],
      ni = all_risk$incidence[j]
    )
   
    all_results = bind_rows(all_results,
    data.frame(
      d = d,
      c = c, 
      age = all_risk$age_group[j],
      ar = i,
      name = "all"
    )
  )

  }
}

trial_results = data.frame(
  d = as.numeric(),
  c = as.numeric(),
  age = as.character(),
  ar = as.integer()
)

for(j in 1:2) {
  for (i in 1:length(AR)) {
    
    d = calculate_prob_freq(
      xi = all_risk$deaths[j] * VE * AR[i],
      xv = all_risk$trial_deaths[j],
      nv = all_risk$trial_pop[j],
      ni = all_risk$incidence[j]
    )
    
    c = calculate_prob_freq(
      xi = all_risk$cases[j] * VE * AR[i],
      xv = all_risk$trial_sae[j],
      nv = all_risk$trial_pop[j],
      ni = all_risk$incidence[j]
    )
    
    trial_results = bind_rows(trial_results,
                        data.frame(
                          d = d, 
                          c = c, 
                          age = all_risk$age_group[j],
                          ar = i,
                          name = "trial"
                        )
    )
    
  }
}

results = bind_rows(all_results, trial_results)


out = results %>% 
  pivot_longer(d:c, names_to = "outcome") %>% 
  mutate(scenario = factor(ar, levels = c(2,3,4,1),
                           labels = c("Small outbreak",
                                      "Large outbreak",
                                      "Endemic",
                                      "Traveller"))) %>%  
  mutate(scenario = fct_rev(scenario)) %>% 
  ggplot(aes(x = 1- value, y = age)) +
  geom_point(aes(color = scenario, shape = name), size = 2) +
  xlim(0,1)  +
  my_theme +
  scale_color_manual(values= cols) +
  scale_shape_manual(values= c(16, 5)) +
  labs(y = "Age group", x = "") +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))  +
  geom_vline(xintercept = 0.95, linetype = 2) +
  facet_wrap(~ outcome)


ggsave(
  out,
  filename = "output/prop_test_figure.jpeg",
  units = "cm",
  height = 6,
  width = 14
)
