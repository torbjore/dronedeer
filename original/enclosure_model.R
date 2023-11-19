# MODEL FOR FENCED DEER
library(nimble)
#setwd("C:/Users/Bruker/Dropbox (UiO)/Deer/Analysis/Model")

# Defining the model
DoubleObsMultisiteCode_fence <- nimbleCode({
  # Model
  for(s in 1:N_surv){
    Psum[s] <- 1-(1-p[s])*(1-p[s])
    pi[1,s] <- p[s]*(1-p[s])/Psum[s]
    pi[2,s] <- (1-p[s])*p[s]/Psum[s]
    pi[3,s] <- p[s]*p[s]/Psum[s]
    mu[s] <- log(mean_lambda) - 0.5*sigma[s]*sigma[s]
    for(i in 1:N_sites[s]){
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      lambda[s,i] <- exp(mu[s] + epsilon[s,i])
      epsilon[s,i] ~ dnorm(0, sd = sigma[s])
      # Observation model:
      Y[s,i] ~ dbin(Psum[s], N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3,s], Y[s,i])
    }
  }
  
  # Priors
  for(s in 1:N_surv){
    sigma[s] ~ dunif(0, 5) # hist(rgamma(1000000, 0.1, 0.1))
    p[s] ~ dunif(0, 1) 
  }
  #mean_lambda ~ dunif(lamblow, lambupp)
  mean_lambda <- 117/5 # Deterministic! Not stochastic prior!
  
  # Derived parameters
  for(s in 1:N_surv){
    median_lambda[s] <- exp(mu[s])
    # For check (should be the same): mean_lambda_site[s] <- exp(mu[s] + 0.5*sigma[s]^2)
  }
  
  # logit_p1 <- log(p1/(1-p1))
  # logit_p2 <- log(p2/(1-p2))
})

# Loading the data
load("../../Analysis/Make_CountData/CountData_Fence.RData")

UseData_fence = CountData_Fence[!is.na(CountData_Fence$Area),]
UseData_fence$Survey = factor(UseData_fence$Survey, levels = unique(UseData_fence$Survey)) # sorted as in data

N_sites = table(UseData_fence$Survey)
UseData_fence$Site = unlist(lapply(N_sites, function(i) 1:i))

Counts_fence = rbind(
  cbind(UseData_fence[,1:4], Site = UseData_fence$Site, Cat = "only_ab", Count = UseData_fence$n_only_ab),
  cbind(UseData_fence[,1:4], Site = UseData_fence$Site, Cat = "only_jb", Count = UseData_fence$n_only_jb),
  cbind(UseData_fence[,1:4], Site = UseData_fence$Site, Cat = "both",    Count = UseData_fence$n_both)
)
Counts_fence$Cat = factor(Counts_fence$Cat, levels = c("only_ab", "only_jb", "both"))

y = tapply(Counts_fence$Count, list(Counts_fence$Survey, Counts_fence$Site, Counts_fence$Cat), sum)
Y = apply(y, c(1,2), sum)
area = tapply(UseData_fence$Area, list(UseData_fence$Survey, UseData_fence$Site), sum)/10000 # Area in hectare


# Function for initial values
colsumy = apply(y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
phat = (p1hat+p2hat)/2
Nhat = apply(Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
lambdahat = (Nhat/N_sites)/apply(area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)

Inits = function(){
  N = round(Y/(1-(1-p1hat)*(1-p2hat)), 0)
  N[is.na(N)] = 0
  sigma = runif(2, 0, 0.1)
  list(
    #mean_lambda = mean(lambdahat)*runif(length(mean(lambdahat)), 0.9, 1.1),
    epsilon = matrix(rnorm(nrow(Y)*ncol(Y), 0, sigma) , nrow = nrow(Y), ncol = ncol(Y)),
    p = exp(log(phat)*runif(2, 0.9, 1.1)),
    N = N,
    sigma = sigma
  )
}

# Setting up the MCMC
DoubleObsMultisiteModel <- nimbleModel(
  DoubleObsMultisiteCode_fence,
  constants = list(N_surv = length(N_sites), N_sites = N_sites), # 0.1 to 10 times point estimate
  data = list(y=y, Y=Y, area = area),
  inits = Inits()
)

CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel)
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("median_lambda", "p", "mu", "sigma"))
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

posterior_samples_fence_p <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=500000,
  nburnin=100000,
  nchain=3,
  thin = 2,
  inits = Inits,
  samplesAsCodaMCMC = TRUE)

library(coda)
plot(posterior_samples_fence_p)

summary(posterior_samples_fence_p)

#save(posterior_samples_fence_p, file = "posterior_samples_fence_p.RData")