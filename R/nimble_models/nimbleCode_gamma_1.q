# Model with gamma distributed lambda and first order covariate effect

# DEFINING THE MODEL
nimbleCode_DOMM_gamma_1 <- nimbleCode({ # Double observer multi-site model (DOMM)

  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      
      # Random p's
      logit_p1[s,i] ~ dnorm(eta1, sd = sigma_p)
      logit_p2[s,i] ~ dnorm(eta2, sd = sigma_p)
      p1[s,i] <- 1/(1+exp(-logit_p1[s,i]))
      p2[s,i] <- 1/(1+exp(-logit_p2[s,i]))
      
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
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta   
      
    }
  }

  # Priors
  for(k in 1:N_sam){
    mu0[k] ~ dnorm(0, sd = 10)
  }

  beta ~ dnorm(0, sd=2)
  
  sigma ~ dgamma(0.01, 0.01)
  shape <- 1/(exp(sigma^2) - 1)
  # Seen with the dunif() prior: If sigma is allowed to be too small, chains sometimes hover around smaller values
  # before jumping to the higher values (but not the other way it looks).
  # I interpret this to mean that there is a local minimum with low sigma.
  
  eta1 ~ dnorm(prior_mean_eta, sd = prior_sd_eta) 
  eta2 ~ dnorm(prior_mean_eta, sd = prior_sd_eta)
  sigma_p ~ T(dgamma(1, 0.05), 0, 0.59)
  
  # Derived parameters
  for(s in 1:N_surv){
    N_tot[s] <- sum(N[s, 1:N_sites[s]])
    sum_area[s] <- sum(area[s, 1:N_sites[s]])
  }
  for(k in 1:N_sam){
    mean_ds[k] <- sum(N_tot[S[1,k]:S[2,k]])/sum(sum_area[S[1,k]:S[2,k]])
  }
  
  # For posterior predictive checks
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      New_N[s,i] ~ dpois(New_lambda[s,i]*area[s,i])
      New_lambda[s,i] ~ dgamma(shape, rate[s,i])
      New_Y[s,i] ~ dbin(Psum[s,i], New_N[s,i])
      New_y[s,i,1:3] ~ dmulti(pi[s,i,1:3], Y[s,i])
      E_Y[s,i] <- Psum[s,i]*area[s,i]*shape/rate[s,i] + 0.0001
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
