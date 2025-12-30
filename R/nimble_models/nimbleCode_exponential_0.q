# Model with exponential lambda and no covariate effect

# DEFINING THE MODEL
nimbleCode_DOMM_exponential_0 <- nimbleCode({
  # Model
  
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
      lambda[s,i] ~ dexp(lamexp[s,i])
      lamexp[s,i] <- 1/exp(mu[s,i])
      mu[s,i] <- mu0[sam[s]]
    }
  }
  
  # Priors
  for(k in 1:N_sam){
    mu0[k] ~ dnorm(0, sd = 10)
  }
  
  eta1 ~ dnorm(prior_mean_eta, sd = prior_sd_eta) 
  eta2 ~ dnorm(prior_mean_eta, sd = prior_sd_eta)
  sigma_p ~ T(dgamma(1, 0.05), 0, 0.59)
  
  # Derived:  mean_ds
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
      New_lambda[s,i] ~ dexp(lamexp[s,i])
      New_N[s,i] ~ dpois(New_lambda[s,i]*area[s,i])
      New_Y[s,i] ~ dbin(Psum[s,i], New_N[s,i])
      New_y[s,i,1:3] ~ dmulti(pi[s,i,1:3], Y[s,i])
      E_Y[s,i] <- Psum[s,i]*area[s,i]/exp(mu[s,i]) + 0.0001
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

