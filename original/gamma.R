setwd("C:/Users/Bruker/Dropbox (UiO)/Deer/Analysis/Model")

# R-packages used: 
library(nimble)
library(coda)

# LOADING MODIFIED PRIORS AND DATA
load(file = "Manuscript_models/mu_logit_p.RData")
load(file = "Manuscript_models/sigma_logit_p.RData")

load(file = "Manuscript_models/UseData.RData")
N_sites = table(UseData$Survey)

load(file = "Manuscript_models/Counts.RData")

# Making y, Y, area and mean_field_dist 
y = tapply(Counts$Count, list(Counts$Survey, Counts$Site, Counts$Cat), sum) # 3D array rekke = survey, colonne = site (bilde), lag = Cat
Y = apply(y, c(1,2), sum) # Matrix with survey x no. deer
area = tapply(UseData$Area, list(UseData$Survey, UseData$Site), sum)/10000 # Area in hectare
mean_field_dist = tapply(UseData$mean_field_dist, list(UseData$Survey, UseData$Site), sum)

# Making standardize function
standardize = function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

# Making index for study area and month
Counts$StudyArea_Month = paste(Counts$Study_area, Counts$Month, sep="_")
tab = table(Counts$Survey,Counts$StudyArea_Month)
sam = apply(tab,1, function(i) which(i != 0))


# DEFINING THE MODEL
DoubleObsMultisiteCode <- nimbleCode({
  # Model
  Psum <- 1-(1-p1)*(1-p2)
  pi[1] <- p1*(1-p2)/Psum
  pi[2] <- (1-p1)*p2/Psum
  pi[3] <- p1*p2/Psum
  
  for(s in 1:N_surv){
    for(i in 1:N_sites[s]){
      # Process model:
      N[s,i] ~ dpois(lambda[s,i]*area[s,i])
      #lambda[s,i] <- exp(mu[s,i] + epsilon[s,i])
      lambda[s,i] ~ dgamma(shape[s,i], rate[sam[s]])
      shape[s,i] <- exp(mu[s,i])*rate[sam[s]]
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta
      #epsilon[s,i] ~ dnorm(0, sd = sigma)
      # Observation model:
      Y[s,i] ~ dbin(Psum, N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3], Y[s,i])
    }
  }
  
  # Priors
  sigma ~ dgamma(0.1, 0.1) 
  
  for(k in 1:N_sam){
    mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
    mean_lambda[k] <- exp(mu0[k])
    mean_lnorm[k] <- log(mean_lambda[k]) - 0.5*sigma*sigma
    var_lognormal[k] <- (exp(sigma*sigma) - 1)*exp(2*mean_lnorm[k] + sigma*sigma)  
    rate[k] <- mean_lambda[k]/var_lognormal[k]
  }
  
  beta ~ dnorm(0, sd=2)
  
  logit_p1 ~ dnorm(mu_logit_p, sd = sigma_logit_p) 
  logit_p2 ~ dnorm(mu_logit_p, sd = sigma_logit_p)
  p1 <- exp(logit_p1)/(1+exp(logit_p1))
  p2 <- exp(logit_p2)/(1+exp(logit_p2))
})


# CREATING FUNCTION FOR INITIAL VALUES
colsumy = apply(y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
Nhat = apply(Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat))
lambdahat_surv = (Nhat/N_sites)/apply(area, 1, mean, na.rm=TRUE) + 0.01 
lambdahat = tapply(lambdahat_surv, list(sam), mean)

Inits = function(){
  N = round(Y/(1-(1-p1hat)*(1-p2hat)), 0)
  N[is.na(N)] = 0
  sigma = runif(1, 1, 2)
  p1 = exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 = exp(log(p2hat)*runif(1, 0.9, 1.1))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    logit_p1 = rnorm(mu_logit_p, sigma_logit_p/5),
    logit_p2 = rnorm(mu_logit_p, sigma_logit_p/5),
    N = N,
    sigma = sigma,
    beta = runif(1, -0.5, 0.5)
  )
}


# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  DoubleObsMultisiteCode,
  constants = list(lamblow = 0.1*lambdahat,  # 0.1 to 10 times point estimate
                   lambupp = 10*lambdahat,
                   N_surv = length(N_sites),
                   N_sites = N_sites,
                   mu_logit_p = mu_logit_p,
                   sigma_logit_p = sigma_logit_p,
                   sam = sam,
                   N_sam = length(unique(sam)),
                   area = area,
                   x = standardize(mean_field_dist)
  ),
  data = list(y=y, Y=Y),
  inits = Inits()
)

CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel)
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("mean_lambda" , "p1", "p2", "mu0", "sigma", "beta"), enableWAIC = TRUE)
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

posterior_gamma <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=500000,
  nburnin=100000,
  nchain=3,
  thin=2,
  inits = Inits,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)

plot(posterior_gamma$samples)
summary(posterior_gamma$samples)
posterior_gamma$WAIC

save(posterior_gamma, file = "Manuscript_models/posterior_gamma.RData")