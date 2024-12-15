# Double observer multi-site model with lognormal lambda and first order covariate effect

# DEFINING THE MODEL
nimbleCode_DOMM_lognormal_2 <- nimbleCode({
  # Model
  Psum <- 1-(1-p1)*(1-p2)
  pi[1] <- p1*(1-p2)/Psum
  pi[2] <- (1-p1)*p2/Psum
  pi[3] <- p1*p2/Psum
  
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      lambda[s,i] <- exp(mu[s,i] + epsilon[s,i])
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta[1] + x[s,i]^2*beta[2]
      epsilon[s,i] ~ dnorm(0, sd = sigma)
      # Observation model:
      Y[s,i] ~ dbin(Psum, N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3], Y[s,i])
    }
  }
  
  # Priors
  for(k in 1:N_sam){
    mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
  }
  
  for(i in 1:2){  
    beta[i] ~ dnorm(0, sd=2) # assume that x is standardized (x_st = (x-mean(x))/sd(x))
  }
  
  logit_p1 ~ dnorm(mu_logit_p, sd = sigma_logit_p) 
  logit_p2 ~ dnorm(mu_logit_p, sd = sigma_logit_p)
  p1 <- exp(logit_p1)/(1+exp(logit_p1))
  p2 <- exp(logit_p2)/(1+exp(logit_p2))
  
  sigma ~ dunif(0.2, 3)

  # Derived: Mean and median densities at mean x
  for(k in 1:N_sam){
    median_lambda[k] <- exp(mu0[k])
    mean_lambda[k] <- exp(mu0[k] + 0.5*sigma^2)
  }
  
  # For posterior predictive checks
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      New_Y[s,i] ~ dbin(Psum, N[s,i])
      New_y[s,i,1:3] ~ dmulti(pi[1:3], Y[s,i])
      E_Y[s,i] <- Psum*N[s,i] + 0.0001
      # Discrepancy contributions
      DiscC_Y[s,i] <- (Y[s,i] - E_Y[s,i])^2 / E_Y[s,i]
      DiscC_New_Y[s,i] <- (New_Y[s,i] - E_Y[s,i])^2 / E_Y[s,i]
      for(j in 1:3){
        E_y[s,i,j] <- pi[j]*Y[s,i] + 0.0001
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
