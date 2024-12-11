# MODEL FOR FENCED DEER
# Modified from original to estimate N instead of density

library(nimble)

# Defining the model
DoubleObsMultisiteCode_fence <- nimbleCode({
  # Model
  for(s in 1:N_surv){
    Psum[s] <- 1-(1-p1[s])*(1-p2[s])
    pi[1,s] <- p1[s]*(1-p2[s])/Psum[s]
    pi[2,s] <- (1-p1[s])*p2[s]/Psum[s]
    pi[3,s] <- p1[s]*p2[s]/Psum[s]
    #mu[s] <- log(mean_lambda) - 0.5*sigma[s]*sigma[s]
    for(i in 1:N_sites[s]){
      # Process model:
      #N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      #lambda[s,i] <- exp(mu[s] + epsilon[s,i])
      #epsilon[s,i] ~ dnorm(0, sd = sigma[s])
      # Observation model:
      Y[s,i] ~ dbin(Psum[s], N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3,s], Y[s,i])
    }
  }
  
  # Priors
  for(s in 1:N_surv){
    #sigma[s] ~ dunif(0, 5) # hist(rgamma(1000000, 0.1, 0.1))
    p1[s] ~ dunif(0, 1)
    p2[s] ~ dunif(0, 1)
    for(i in 1:N_sites[s]){
      N[s,i] ~ dpois(lambda_N[s,i]) # dflat() # dunif(0, 100)
      lambda_N[s,i] ~ dunif(0, 50)
    }
  }
  #lambda_N ~ dunif(0, 35)
  #mean_lambda ~ dunif(lamblow, lambupp)
  #mean_lambda <- 117/5 # Deterministic! Not stochastic prior!
  
  # Derived parameters
  for(s in 1:N_surv){
    N_tot[s] <- sum(N[s, 1:N_sites[s]])
    Dens[s] <- N_tot[s]/sum(area[s, 1:N_sites[s]])
  }
  Dens_tot <- sum(N_tot[1:2])/(sum(area[1, 1:N_sites[1]]) + sum(area[2, 1:N_sites[2]]))
  
  
  # for(s in 1:N_surv){
  #   median_lambda[s] <- exp(mu[s])
  #   # For check (should be the same): mean_lambda_site[s] <- exp(mu[s] + 0.5*sigma[s]^2)
  # }
  
  # logit_p1 <- log(p1/(1-p1))
  # logit_p2 <- log(p2/(1-p2))
})