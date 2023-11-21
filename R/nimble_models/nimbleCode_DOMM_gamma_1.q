# Double observer multi-site model with lognormal lambda and first order covariate effect

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
      shape[s,i] <- exp(mu[s,i])*rate
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta
      # Observation model:
      Y[s,i] ~ dbin(Psum, N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3], Y[s,i])
    }
  }
  
  # Priors
  #rate ~ dunif(0,10)
  rate ~ dgamma(0.1, 0.1)
  
  # sigma ~ dgamma(0.1, 0.1) 
  # 
  for(k in 1:N_sam){
  mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
  mean_lambda[k] <- exp(mu0[k])
  #   mean_lnorm[k] <- log(mean_lambda[k]) - 0.5*sigma*sigma
  #   var_lognormal[k] <- (exp(sigma*sigma) - 1)*exp(2*mean_lnorm[k] + sigma*sigma)  
  #   rate[k] <- mean_lambda[k]/var_lognormal[k]
  }
  # 
  beta ~ dnorm(0, sd=2)
  
  logit_p1 ~ dnorm(mu_logit_p, sd = sigma_logit_p) 
  logit_p2 ~ dnorm(mu_logit_p, sd = sigma_logit_p)
  p1 <- exp(logit_p1)/(1+exp(logit_p1))
  p2 <- exp(logit_p2)/(1+exp(logit_p2))
})
