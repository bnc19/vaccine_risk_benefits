library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

my_theme = theme_bw() +
  theme(
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.title = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )

cols = c(  "#FD8D3C", 
           "#F768A1",
           "#C51B8A",
           "#7A0177" )

cols2 = c("skyblue", "#FDB863", "#49006A")

# Generate toy data for U-shaped infection risk
ages = seq(0, 100, by = 1)

# U-shaped infection risk (high in neonates and elderly)
infection_risk = 0.02 + 0.3 * exp(-ages/5) + 0.2 * exp((ages - 100)/20)

# Two vaccine risk scenarios
vaccine_risk_constant = rep(0.006, length(ages))  # Independent of age
vaccine_risk_age_dependent = 0.001 + 0.08 * (ages/100)^2  # Increases with age

# Create data frame for infection and vaccine risks
risk_data = data.frame(
  age = ages,
  infection_risk = infection_risk,
  vaccine_risk_constant = vaccine_risk_constant,
  vaccine_risk_age_dependent = vaccine_risk_age_dependent
)

# Plot 1: U-shaped infection risk

p1 = risk_data %>% 
  pivot_longer(cols = -age, names_to = "risk_type", values_to = "risk") %>%
  mutate(risk_type = recode(risk_type,
                            "infection_risk" = "Age-dependent infection risk",
                            "vaccine_risk_constant" = "Constant vaccine risk",
                            "vaccine_risk_age_dependent" = "Age-dependent vaccine risk")) %>% 
  ggplot(aes(x = age, y = risk)) +
  geom_line(aes(color = risk_type), linewidth = 0.7) +
  my_theme + 
  theme(legend.position = c(0.61,0.88),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm")) +
  labs(x = "Age (years)",
       y = "Probability of adverse outcome") +
  ylim(0,0.35)+
  scale_color_manual(values=cols2)


print(p1)

# Plot 2: Risk-benefit analysis with baseline efficacy
VE = 0.95
AR = 0.3
pop = 10000

# Corrected the risk calculations with proper parentheses
risk_nv =  AR * risk_data$infection_risk * pop
risk_v_c = AR * risk_data$infection_risk * (1-VE) + risk_data$vaccine_risk_constant * pop
risk_v_a = AR * risk_data$infection_risk * (1-VE) + risk_data$vaccine_risk_age_dependent * pop
benefit_c = risk_nv - risk_v_c
benefit_a = risk_nv - risk_v_a 

risk_benefit_data = risk_data %>%
  mutate(
    risk_nv = risk_nv,
    risk_v_c = risk_v_c,
    risk_v_a = risk_v_a,
    benefit_constant = benefit_c,
    benefit_age_dependent = benefit_a,
  )

risk_benefit_long = risk_benefit_data %>%
  select(age, benefit_constant, benefit_age_dependent) %>%
  pivot_longer(cols = -age, names_to = "vaccine_type", values_to = "benefit") %>%
  mutate(vaccine_type = recode(vaccine_type,
                               "benefit_constant" = "Constant vaccine risk",
                               "benefit_age_dependent" = "Age-dependent vaccine risk"))

p2 = ggplot(risk_benefit_long, aes(x = age, y = benefit)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(aes(color = vaccine_type), linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  my_theme + 
  labs(x = "Age (years)",
       y = "Adverse outcomes\naverted (per 10,000)") +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  scale_color_manual(values = c("#FDB863", "#49006A")) +
  scale_y_continuous(
    limits = c(min(risk_benefit_long$benefit, na.rm = TRUE), 
               max(risk_benefit_long$benefit, na.rm = TRUE)),
    breaks = c(-200,0,200,400,600,800,1000)
  )

print(p2)

# Plot 4: Risk-benefit varies with efficacy
efficacy_levels = seq(0, 1, by = 0.1)
age_points = c(1, 20, 60, 80) 

efficacy_data = expand.grid(efficacy = efficacy_levels, age = age_points)

# Calculate risk-benefit for different efficacies
efficacy_data = efficacy_data %>%
  left_join(risk_data, by = "age") %>%
  mutate(
    benefit_constant = (AR * infection_risk * pop) - ((AR * infection_risk * (1-efficacy) + vaccine_risk_constant) *  pop)
  )

# Plot for constant vaccine risk
p3 = ggplot(efficacy_data, aes(x = efficacy*100, y = benefit_constant, color = factor(age))) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  my_theme +
  labs(x = "Vaccine efficacy (100%)",
       color = "Age (years)",
       y = "Adverse outcomes\naverted (per 10,000)") +
  theme(legend.position = c(0.17, 0.75),
        legend.key.size = unit(0.3, "cm"))+
  scale_color_manual(values = cols)+
  scale_y_continuous(
    limits = c(min(efficacy_data$benefit, na.rm = TRUE), 
               max(efficacy_data$benefit, na.rm = TRUE)),
    breaks = c(-60,0,200,400,600, 800)
  )


print(p3)

# Plot 5: Risk-benefit varies with attack rate
attack_rates = seq(0, 1, by = 0.1)

attack_data = expand.grid(attack_rate = attack_rates, age = age_points)

# Calculate risk-benefit for different attack rates
attack_data = attack_data %>%
  left_join(risk_data, by = "age") %>%
  mutate(
    benefit_constant = (attack_rate * infection_risk * pop) - ((attack_rate * infection_risk * (1-0.95) + vaccine_risk_constant) * pop)
  )

# Plot for constant vaccine risk
p4 = ggplot(attack_data, aes(x = attack_rate*100, y = benefit_constant, color = factor(age))) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.5) +
  geom_line(size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  my_theme +
  labs(x = "Attack rate (%)",
       y = "Adverse outcomes\naverted (per 10,000)", 
       color = "Age")  +
  theme(legend.position = "none") +
  scale_color_manual(values = cols) +
  scale_y_continuous(
    limits = c(-200, 
               max(attack_data$benefit, na.rm = TRUE)),
    breaks = c(-200,0,500,1000,1500,2000,2500)
  )

print(p4)


# save all plots 

out = wrap_plots(
  wrap_plots(p1, p2,  ncol = 2),
  wrap_plots(p3, p4 , ncol = 2),
  ncol = 1
) +
  plot_annotation(tag_levels = "a")

out

ggsave(out, 
       filename = "output/example_risk.jpg",
       units = "cm",
       height = 12, width = 14)

