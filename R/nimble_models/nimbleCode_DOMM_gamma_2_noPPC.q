# Double observer multi-site model with gamma distributed lambda and second order covariate effect
# POSTERIOR PREDICTIVE CHECKS REMOVED FOR SPEED

# DEFINING THE MODEL
nimbleCode_DOMM_gamma_2 <- nimbleCode({

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
      lambda[s,i] ~ dgamma(shape, rate[s,i])
      rate[s,i] <- shape/exp(mu[s,i] + 0.5*sigma^2)
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta[1] + x[s,i]^2*beta[2]     
      
    }
  }

  # Priors
  for(k in 1:N_sam){
    # mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
    mu0[k] ~ dnorm(mu0hat[k], sd = 1.175)
  }
  
  for(i in 1:2){  
    beta[i] ~ dnorm(0, sd=2) # assume that x is standardized (x_st = (x-mean(x))/sd(x))
  }
  
  sigma ~ dunif(0.5, 3)
  shape <- 1/(exp(sigma^2) - 1)
  # If sigma is allowed to be too small, chains sometimes hover around smaller values
  # before jumping to the higher values (but not the other way it looks).
  # I interpret this to mean that there is a local minimum with low sigma.
  
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
  # chains entered flat parts of the likelihood surface, we also set a lower
  # bound of 0.1.
  })
