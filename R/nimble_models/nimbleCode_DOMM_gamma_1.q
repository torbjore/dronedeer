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
      lambda[s,i] ~ dgamma(shape, rate[s,i])
      rate[s,i] <- shape/exp(mu[s,i] + 0.5*sigma^2)
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta      
      
      # Observation model:
      Y[s,i] ~ dbin(Psum, N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3], Y[s,i])
    }
  }

  # Negative binomial variant
  # P ~ dunif(0,1)
  # rate <- P/(1-P)
    
  # # Priors
  # invrate ~ dunif(0, 5)# 10)
  # rate <- 1/invrate
  # # NB! It may be that we need an upper bound on rate if there is a local peak
  # with small variance. Could maybe also try a uniform prior for 1/rate.

  for(k in 1:N_sam){
    #exp_mu0[k] ~ dunif(0, 10) # Denne bør nok settes mer restrektivt, slik som i log-normal modellen
    exp_mu0[k] ~ dunif(lamblow[k], lambupp[k])
    mu0[k] <- log(exp_mu0[k])
  }

  sigma ~ dunif(0.5, 3)
  shape <- 1/(exp(sigma^2) - 1)
  
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
