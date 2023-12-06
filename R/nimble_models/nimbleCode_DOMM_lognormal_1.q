# Double observer multi-site model with lognormal lambda and first order covariate effect

# DEFINING THE MODEL
nimbleCode_DOMM_lognormal_1 <- nimbleCode({
  # Model
  
  # Fixed p's
  # Psum <- 1-(1-p1)*(1-p2)
  # pi[1] <- p1*(1-p2)/Psum
  # pi[2] <- (1-p1)*p2/Psum
  # pi[3] <- p1*p2/Psum
  
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      # Random p's
      logit_p1[s,i] ~ dnorm(mu_p1, sd = sigma_p)
      logit_p2[s,i] ~ dnorm(mu_p2, sd = sigma_p)
      p1[s,i] <- exp(logit_p1[s,i])/(1+exp(logit_p1[s,i]))
      p2[s,i] <- exp(logit_p2[s,i])/(1+exp(logit_p2[s,i]))
      
      Psum[s,i] <- 1-(1-p1[s,i])*(1-p2[s,i])
      pi[s,i,1] <- p1[s,i]*(1-p2[s,i])/Psum[s,i]
      pi[s,i,2] <- (1-p1[s,i])*p2[s,i]/Psum[s,i]
      pi[s,i,3] <- p1[s,i]*p2[s,i]/Psum[s,i]
      
      # Observation model:
      Y[s,i] ~ dbin(Psum[s,i], N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[s, i, 1:3], Y[s,i])
      
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      lambda[s,i] <- exp(mu[s,i] + epsilon[s,i])
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta
      epsilon[s,i] ~ dnorm(0, sd = sigma)
    }
  }
  
  # Priors
  for(k in 1:N_sam){
    mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
  }
  
  beta ~ dnorm(0, sd=2) # assume that x is standardized (x_st = (x-mean(x))/sd(x))
  
  mu_p1 ~ dnorm(prior_mu_logit_p, sd = prior_sigma_logit_p) 
  mu_p2 ~ dnorm(prior_mu_logit_p, sd = prior_sigma_logit_p)
  sigma_p ~ dunif(0.1, 0.59)
  # It is not reasonable that sigma_p could be very high because observers did an
  # honest attempt to locate deer at all sites, and we we don't believe that a
  # very high proportion of sites with zero counts actually had deer present (it's
  # more likely that most of these sites did not have deer in them). Assuming
  # that the 97.% percentile site had a maximum of 10 times higher detection
  # odds as the 2.5% percentile site, we used a uniform prior for sigma_p with 
  # an upper bound of log(10)/(1.96*2) = 0.59. To avoid poor mixing when the 
  # chains entered flatparts of the likelihood surface, we also set a lower
  # bound of 0.1.
  
  sigma ~ dunif(0.5, 3) # WAIC 122.1556 (before random effects on p)
# sigma ~ dunif(0, 0.5) # WAIC 131.1987 (before random effects on p), chains hit the roof, and GOF is worse (Bayesian p-value based on Pearson statistic for Y 0.2989667)
  # If allowed to be too small, chains sometimes hover around smaller values
  # before jumping to the higher values (but not the other way it looks).
  # I interpret this to mean that there is a local minimum with low sigma.
  
  # Derived: Mean and median densities at mean x
  for(k in 1:N_sam){
    median_lambda[k] <- exp(mu0[k])
    mean_lambda[k] <- exp(mu0[k] + 0.5*sigma^2)
  }
  
  # For posterior predictive checks
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      New_Y[s,i] ~ dbin(Psum[s,i], N[s,i])
      New_y[s,i,1:3] ~ dmulti(pi[s,i,1:3], Y[s,i])
      E_Y[s,i] <- Psum[s,i]*N[s,i] + 0.0001
      # Discrepancy contributions
      DiscC_Y[s,i] <- (Y[s,i] - E_Y[s,i])^2 / E_Y[s,i]
      DiscC_New_Y[s,i] <- (New_Y[s,i] - E_Y[s,i])^2 / E_Y[s,i]
      for(j in 1:3){
        E_y[s,i,j] <- pi[s,i,j]*Y[s,i] + 0.00001
        DiscC_y_v[s,i, j] <- (y[s,i,j] - E_y[s,i,j])^2 / E_y[s,i,j]
        DiscC_New_y_v[s,i, j] <- (New_y[s,i,j] - E_y[s,i,j])^2 / E_y[s,i,j]
      }
      DiscC_y[s,i] <- sum(DiscC_y_v[s,i, 1:3])
      DiscC_New_y[s,i] <- sum(DiscC_New_y_v[s,i, 1:3])
    }
    # Sums per survey
    Disc_Y_s[s] <- sum(DiscC_Y[s, 1:N_sites[s]])
    Disc_New_Y_s[s] <- sum(DiscC_New_Y[s, 1:N_sites[s]])
    Disc_y_s[s] <- sum(DiscC_y[s, 1:N_sites[s]])
    Disc_New_y_s[s] <- sum(DiscC_New_y[s, 1:N_sites[s]])
  }
  # Discrepancy statistics (Pearson chi-sq, but note many small expectations)
  Disc_Y <- sum(Disc_Y_s[1:N_surv])
  Disc_New_Y <- sum(Disc_New_Y_s[1:N_surv])
  Disc_y <- sum(Disc_y_s[1:N_surv])
  Disc_New_y <- sum(Disc_New_y_s[1:N_surv])
})
