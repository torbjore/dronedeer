# MODEL FOR FENCED DEER

library(nimble)

# Defining the model
DoubleObsMultisiteCode_fence <- nimbleCode({
  # Model
  
  for(s in 1:N_surv){
    logit_p1[s] ~ dnorm(mu_p, sd = sigma_p)
    logit_p2[s] ~ dnorm(mu_p, sd = sigma_p)
    p1[s] <- exp(logit_p1[s])/(1+exp(logit_p1[s]))
    p2[s] <- exp(logit_p2[s])/(1+exp(logit_p2[s]))
    Psum[s] <- 1-(1-p1[s])*(1-p2[s])
    pi[s,1] <- p1[s]*(1-p2[s])/Psum[s]
    pi[s,2] <- (1-p1[s])*p2[s]/Psum[s]
    pi[s,3] <- p1[s]*p2[s]/Psum[s]
    mu[s] <- log(mean_lambda) - 0.5*sigma[s]*sigma[s]
    for(i in 1:N_sites[s]){
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      lambda[s,i] <- exp(mu[s] + epsilon[s,i])
      epsilon[s,i] ~ dnorm(0, sd = sigma[s])
      # Observation model:
      Y[s,i] ~ dbin(Psum[s], N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[s,1:3], Y[s,i])
    }
  }
  
  # Priors
  for(s in 1:N_surv){
    sigma[s] ~ dunif(0, 5)
    #p[s] ~ dunif(0, 1)
    # for(i in 1:N_sites[s]){
    #   N[s,i] ~ dpois(lambda_N[s,i]) # dflat() # dunif(0, 100)
    #   lambda_N[s,i] ~ dunif(0, 50)
    # }
  }
  #lambda_N ~ dunif(0, 35)
  mean_lambda <- 117/5 # Deterministic! Not stochastic prior!
  
  mu_p ~ dnorm(0, sd = 2) 
  sigma_p ~ dunif(0.1, 0.59)
  
  # Derived parameters
  for(s in 1:N_surv){
    N_tot[s] <- sum(N[s, 1:N_sites[s]])
    Dens[s] <- N_tot[s]/sum(area[s, 1:N_sites[s]])
  }
  Dens_tot <- sum(N_tot[1:2])/(sum(area[1, 1:N_sites[1]]) + sum(area[2, 1:N_sites[2]]))
  
  # For posterior predictive checks
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      New_Y[s,i] ~ dbin(Psum[s], N[s,i])
      New_y[s,i,1:3] ~ dmulti(pi[s,1:3], Y[s,i])
      E_Y[s,i] <- Psum[s]*N[s,i] + 0.0001
      # Discrepancy contributions
      DiscC_Y[s,i] <- (Y[s,i] - E_Y[s,i])^2 / E_Y[s,i]
      DiscC_New_Y[s,i] <- (New_Y[s,i] - E_Y[s,i])^2 / E_Y[s,i]
      for(j in 1:3){
        E_y[s,i,j] <- pi[s,j]*Y[s,i] + 0.00001
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
  
  # for(s in 1:N_surv){
  #   median_lambda[s] <- exp(mu[s])
  #   # For check (should be the same): mean_lambda_site[s] <- exp(mu[s] + 0.5*sigma[s]^2)
  # }
  
  # logit_p1 <- log(p1/(1-p1))
  # logit_p2 <- log(p2/(1-p2))
})