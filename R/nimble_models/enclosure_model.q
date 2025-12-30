# MODEL FOR FENCED DEER

library(nimble)

# Defining the model
DoubleObsMultisiteCode_fence <- nimbleCode({
  # Model
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      logit_p1[s,i] ~ dnorm(eta[s], sd = sigma_p)
      logit_p2[s,i] ~ dnorm(eta[s], sd = sigma_p)
      p1[s,i] <- 1/(1+exp(-logit_p1[s,i]))
      p2[s,i] <- 1/(1+exp(-logit_p2[s,i]))
      Psum[s,i] <- 1-(1-p1[s,i])*(1-p2[s,i])
      pi[s,i,1] <- p1[s,i]*(1-p2[s,i])/Psum[s,i]
      pi[s,i,2] <- (1-p1[s,i])*p2[s,i]/Psum[s,i]
      pi[s,i,3] <- p1[s,i]*p2[s,i]/Psum[s,i]
      Y[s,i] ~ dbin(Psum[s,i], N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[s,i,1:3], Y[s,i])
    }
  }
  
  # Priors
  for(s in 1:N_surv){
    eta[s] ~ dnorm(0, sd = 2) 
    for(i in 1:N_sites[s]){
      N[s,i] ~ dpois(lambda_N[s,i]*area[s,i]) # hierarchical prior
      lambda_N[s,i] ~ dunif(0, 50)
    }
  }
  sigma_p ~ T(dgamma(1, 0.05), 0, 0.59)
  
  # Derived parameters
  for(s in 1:N_surv){
    N_sum[s] <- sum(N[s, 1:N_sites[s]]) # This is the total in the sites surveyed, not the entire enclosure
    #Dens[s] <- N_sum[s]/sum_area[s]
    N_tot[s] <- N_sum[s]/(sum_area[s]/5) # Enclosure is 5 ha
    var_N_tot[s] <- (N_sites[s]/(sum_area[s]/5) - N_sites[s])*(sd(N[s, 1:N_sites[s]])^2)/(sum_area[s]/5)
    # lower and upper CI based on var_N_tot[s] and assuming normal distribution on a log-scale (central limit theorem)
    N_lower[s] <- exp(log(N_tot[s]) - 2*(N_tot[s]^(-1))*sqrt(var_N_tot[s]))
    N_upper[s] <- exp(log(N_tot[s]) + 2*(N_tot[s]^(-1))*sqrt(var_N_tot[s]))
  }
 
  # For posterior predictive checks
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      New_Y[s,i] ~ dbin(Psum[s,i], N[s,i])
      New_y[s,i,1:3] ~ dmulti(pi[s,i,1:3], Y[s,i])
      E_Y[s,i] <- Psum[s,i]*lambda_N[s,i]*area[s,i] + 0.0001
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