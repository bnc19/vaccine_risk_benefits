## A Framework for Risk-Benefit Analysis of Vaccines Approved Through Accelerated Pathways

This repository contains the data and code needed to reproduce the results presented in the manuscript: A Framework for Risk-Benefit Analysis of Vaccines Approved Through Accelerated Pathways.


Scripts can be run in order of 1-11 to reproduce the papers findings and figures. The purpose of each script are as follows: 

* "1.example_risk.R": Script to create Figure S1, demonstrating the factors affecting the risks and benefits of vaccination. 
* "2.risk_benefit_framework.R": Script to run framework to calculate the number of cases averted given different risks of severe outcomes following infection or vaccine, epidemiological scenarios and vaccine efficacies. File outputs Figure 1. 
* "3.general_prob_benefits.R": Script to calculate the probability that a vaccine benefits outweigh the risks of vaccine SAEs under a Bayesian framework, given different probabilities of severe outcomes following infection, probabilities of adverse vaccine events and either 1,000 or 10,000 individuals vaccinated (cohort size), epidemiological scenarios and vaccine efficacies. File outputs Figure S2. 
* "4.estimating_sample_size.R": Script to calculate sample size needed to achieve 80% power to detect at least one death as a function of attack rate and probability of death from infection. File outputs Figure S5. 
* "5.chik_data.R": Script to plot the available data on IXCHIQ doses, SAEs, and risk of severe outcomes given chikungunya infection. File outputs Figure S3. 
* "6.chik_risk_benefits.R": Script to calculate the expected number of medically attended chikungunya cases / SAEs and deaths per 10,000 individuals in each age group in the presence or absence of IXCHIQ vaccination across different epidemiological scenarios and age groups. File outputs Figure 3A,B,D,E.
* "7.chik_risk_benefits_sensitivity.R": Script to run sensitivity analyses on the impact of VE and AR on expected number of medically attended chikungunya cases / SAEs and deaths averted by IXCHIQ vaccination. File outputs Figures S4 and S7. 
* "8.chik_prob_benefits.R": Script to calculate the probability that IXCHIQs benefits are greater than its risks. File outputs Figures 4 and S6, S7, S8. 
* "9.calculate_dalys.R": Script to calculate DALYs averted by IXCHIQ vaccination. File outputs Figure 3C/F.
* "10.threshold_vaccine_risk.R": # Script to calculate the maximum acceptable vaccine risk for different diseases, attack rates, vaccine efficacies. File outputs Figure 2.
* "11.calculate_dalys.R": Script to calculate DALYs averted by IXCHIQ vaccination if all 3 deaths were vaccine linked. File outputs Figure S9.
* "functions.R" Source functions. 
