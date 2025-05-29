
# Function to extract Beta parameters from binconf results
extract_beta_params = function(x, n, alpha0 = 1, beta0 = 1) {
  # With uniform prior, Beta posterior has parameters:
  alpha = x + alpha0
  beta = n - x + beta0
  return(list(alpha = alpha, beta = beta))
}


simulate_vaccine_benefit = function(VE, AR, all_risk, di, dv, ci, cv, n, pop) {
  
  # Initialize storage for results
  
  summary_results <- tibble()
  
  # Nested loops over i, j, and k
  for (i in 1:length(VE)) {
    for (j in 1:nrow(all_risk)) {
      for (k in 1:length(AR)) {
        
        
        # Generate random samples
        pdi = rbeta(n, shape1 = di$alpha[j], shape2 = di$beta[j])
        pdv = rbeta(n, shape1 = dv$alpha[j], shape2 = dv$beta[j])
        pci = rbeta(n, shape1 = ci$alpha[j], shape2 = ci$beta[j])
        pcv = rbeta(n, shape1 = cv$alpha[j], shape2 = cv$beta[j])
        
        # Calculate threshold and prob benefit
        d_threshold = pdi * VE[i] * AR[k]
        c_threshold = pci * VE[i] * AR[k]
        d_benefit = d_threshold > pdv
        c_benefit = c_threshold > pcv
        
        # Calculate proportion averted 
        d_av = (d_threshold - pdv) * pop
        c_av = (c_threshold - pcv) * pop
        
        # Calculate cases with and without vaccination 
        ndi = pdi * AR[k] * pop 
        nci = pci * AR[k] * pop 
        ncv = nci  * (1-VE[i])  
        ndv = ndi  * (1-VE[i]) 
        tdv = (pdi * AR[k] * (1-VE[i]) + pdv) * pop
        tcv = (pci * AR[k] * (1-VE[i]) + pcv) * pop
        
        VE_value = VE[i]
        AR_value = AR[k]
        age_value = all_risk$age_group[j]
        
        # Store summary results
        summary_results <- bind_rows(summary_results, 
                                     tibble(
                                       ar_scenario = k, 
                                       age_group = age_value,
                                       VE_value = VE_value,
                                       AR_value = AR_value,
                                       d_benefit = mean(d_benefit),
                                       c_benefit = mean(c_benefit),
                                       mean_d_av = mean(d_av),
                                       mean_c_av = mean(c_av),
                                       low_d_av  = quantile(d_av, probs = c(0.025)),
                                       low_c_av  = quantile(c_av, probs = c(0.025)),
                                       high_d_av = quantile(d_av, probs = c(0.975)),
                                       high_c_av = quantile(c_av, probs = c(0.975)),
                                       mean_pd_i = mean(pdi* pop),
                                       mean_pd_v = mean(pdv* pop),
                                       mean_pc_i = mean(pci* pop),
                                       mean_pc_v = mean(pcv* pop),
                                       low_pd_i  = quantile(pdi* pop, probs = c(0.025)),
                                       low_pd_v  = quantile(pdv* pop, probs = c(0.025)),
                                       high_pd_i = quantile(pdi* pop, probs = c(0.975)),
                                       high_pd_v = quantile(pdv* pop, probs = c(0.975)),
                                       low_pc_i  = quantile(pci* pop, probs = c(0.025)),
                                       low_pc_v  = quantile(pcv* pop, probs = c(0.025)),
                                       high_pc_i = quantile(pci* pop, probs = c(0.975)),
                                       high_pc_v = quantile(pcv* pop, probs = c(0.975)),
                                       mean_nd_i = mean(ndi),
                                       mean_nd_v = mean(ndv),
                                       mean_nc_i = mean(nci),
                                       mean_nc_v = mean(ncv),
                                       low_nd_i  = quantile(ndi, probs = c(0.025)),
                                       low_nd_v  = quantile(ndv, probs = c(0.025)),
                                       high_nd_i = quantile(ndi, probs = c(0.975)),
                                       high_nd_v = quantile(ndv, probs = c(0.975)),
                                       low_nc_i  = quantile(nci, probs = c(0.025)),
                                       low_nc_v  = quantile(ncv, probs = c(0.025)),
                                       high_nc_i = quantile(nci, probs = c(0.975)),
                                       high_nc_v = quantile(ncv, probs = c(0.975)),
                                       mean_td_v = mean(tdv),
                                       mean_tc_v = mean(tcv),
                                       low_td_v  = quantile(tdv, probs = c(0.025)),
                                       low_tc_v  = quantile(tcv, probs = c(0.025)),
                                       high_td_v = quantile(tdv, probs = c(0.975)),
                                       high_tc_v = quantile(tcv, probs = c(0.975)),
                                     ))
      }
    }
  }
  
  return(summary_results)
  
}
