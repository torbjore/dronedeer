#setwd("C:/Users/Bruker/Dropbox (UiO)/Deer/Analysis/Model")

# R-packages used: 
library(nimble)
library(coda)

# # LOADING MODIFIED PRIORS AND DATA
# load(file = "Manuscript_models/mu_logit_p.RData")
# load(file = "Manuscript_models/sigma_logit_p.RData")
# 
# load(file = "Manuscript_models/UseData.RData")
# N_sites = table(UseData$Survey)
# 
# load(file = "Manuscript_models/Counts.RData")

# LOADING MODIFIED PRIORS AND DATA
load(file = "mu_logit_p.RData")
load(file = "sigma_logit_p.RData")

load(file = "UseData.RData")
N_sites = table(UseData$Survey)

load(file = "Counts.RData")


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
      lambda[s,i] ~ dexp(lamexp[s,i])
      lamexp[s,i] <- 1/exp(mu[s,i])
      mu[s,i] <- mu0[sam[s]] + x[s,i]*beta
      # Observation model:
      Y[s,i] ~ dbin(Psum, N[s,i])
      y[s, i, 1:3] ~ dmulti(pi[1:3], Y[s,i])
    }
  }
  
  # Priors
  
  for(k in 1:N_sam){
    mu0[k] ~ dunif(log(lamblow[k]), log(lambupp[k]))
    mean_lambda[k] <- exp(mu0[k])
  }
  
  beta ~ dnorm(0, sd=2) # assume that x is standardized (x_standard = (x-mean(x))/sd(x))
  
  logit_p1 ~ dnorm(mu_logit_p, sd = sigma_logit_p) 
  logit_p2 ~ dnorm(mu_logit_p, sd = sigma_logit_p)
  p1 <- exp(logit_p1)/(1+exp(logit_p1))
  p2 <- exp(logit_p2)/(1+exp(logit_p2))
})


# CREATING FUNCTION FOR INITIAL VALUES
colsumy = apply(y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
Nhat = apply(Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
lambdahat_surv = (Nhat/N_sites)/apply(area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)
lambdahat = tapply(lambdahat_surv, list(sam), mean)

Inits = function(){
  N = round(Y/(1-(1-p1hat)*(1-p2hat)), 0)
  N[is.na(N)] = 0
  p1 = exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 = exp(log(p2hat)*runif(1, 0.9, 1.1))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    logit_p1 = rnorm(mu_logit_p, sigma_logit_p/5),
    logit_p2 = rnorm(mu_logit_p, sigma_logit_p/5),
    N = N,
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
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("mean_lambda" , "p1", "p2", "mu0", "beta"), enableWAIC = TRUE)
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

posterior_exponential <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=500000,
  nburnin=100000,
  nchain=3,
  thin=2,
  inits = Inits,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)

plot(posterior_exponentia$samples)
summary(posterior_exponential$samples)
posterior_exponential$WAIC

save(posterior_exponential, file = "Manuscript_models/posterior_exponential.RData")

