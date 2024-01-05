# Double observer multi-site model with gamma distributed lambda and first order covariate effect

# DEFINING THE MODEL
nimbleCode_DOMM_gamma_1 <- nimbleCode({
  # Model
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
      
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      lambda[s,i] ~ dgamma(shape, rate[s,i])
      rate[s,i] <- shape/exp(mu[s,i] + 0.5*sigma^2)
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta      
      
      # Observation model:
      Y[s,i] ~ dbin(Psum[s,i], N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[s, i, 1:3], Y[s,i])
    }
  }

  # for(k in 1:N_sam){
  #   exp_mu0[k] ~ dunif(lamblow[k], lambupp[k])
  #   mu0[k] <- log(exp_mu0[k])
  # }
  
  # Priors
  for(k in 1:N_sam){
    mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
  }

  sigma ~ dunif(0.5, 3)
  shape <- 1/(exp(sigma^2) - 1)
  
  beta ~ dnorm(0, sd=2)
  
  mu_p1 ~ dnorm(prior_mu_logit_p, sd = prior_sigma_logit_p) 
  mu_p2 ~ dnorm(prior_mu_logit_p, sd = prior_sigma_logit_p)
  sigma_p ~ dunif(0.1, 0.59)
  
  # logit_p1 ~ dnorm(mu_logit_p, sd = sigma_logit_p) 
  # logit_p2 ~ dnorm(mu_logit_p, sd = sigma_logit_p)
  # p1 <- exp(logit_p1)/(1+exp(logit_p1))
  # p2 <- exp(logit_p2)/(1+exp(logit_p2))
  
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
