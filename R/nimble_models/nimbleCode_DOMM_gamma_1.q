# Double observer multi-site model with gamma distributed lambda and first order covariate effect

# DEFINING THE MODEL
nimbleCode_DOMM_gamma_1 <- nimbleCode({
  # Model
  Psum <- 1-(1-p1)*(1-p2)
  pi[1] <- p1*(1-p2)/Psum
  pi[2] <- (1-p1)*p2/Psum
  pi[3] <- p1*p2/Psum
  
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      lambda[s,i] ~ dgamma(shape[s,i], rate)
      #shape[s,i] <- exp(mu[s,i])*rate
      #mu[s,i] <- mu0[sam[s]] + x[s,i]*beta
      shape[s,i] <- exp_mu0[sam[s]]*exp(x[s,i]*beta)*rate
      # Observation model:
      Y[s,i] ~ dbin(Psum, N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3], Y[s,i])
    }
  }
  
  # Priors
  rate ~ dunif(0, 10)

  for(k in 1:N_sam){
    exp_mu0[k] ~ dunif(0, 10)
    mu0[k] <- log(exp_mu0[k])
  }

  beta ~ dnorm(0, sd=2)
  
  logit_p1 ~ dnorm(mu_logit_p, sd = sigma_logit_p) 
  logit_p2 ~ dnorm(mu_logit_p, sd = sigma_logit_p)
  p1 <- exp(logit_p1)/(1+exp(logit_p1))
  p2 <- exp(logit_p2)/(1+exp(logit_p2))
  
  # For posterior predictive checks
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      New_Y[s,i] ~ dbin(Psum, N[s,i])
      #New_y[s,i,1:3] ~ dmulti(pi[1:3], Y[s,i])
      E_Y[s,i] <- Psum*N[s,i] + 0.0001
      #E_y[s,i,1:3] <- pi[1:3]*Y[s,i]
      # Discrepancy contributions
      DiscC_Y[s,i] <- pow(Y[s,i] - E_Y[s,i], 2) / E_Y[s,i]
      DiscC_New_Y[s,i] <- pow(New_Y[s,i] - E_Y[s,i], 2) / E_Y[s,i]
    }
    # Sums per survey
    Disc_Y_s[s] <- sum(DiscC_Y[s, 1:N_sites[s]])
    Disc_New_Y_s[s] <- sum(DiscC_New_Y[s, 1:N_sites[s]])
  }
  # Discrepancy statistics (Pearson chi-sq, but note many small expectations)
  Disc_Y <- sum(Disc_Y_s[1:N_surv])
  Disc_New_Y <- sum(Disc_New_Y_s[1:N_surv])
})
