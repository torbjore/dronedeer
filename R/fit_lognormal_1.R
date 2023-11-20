# R-packages used: 
library(nimble)
library(coda)

# Making standardize function (NB! Do not use scale() as this function scales each column of a matrix)
standardize = function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

# KOMMENTAR: legg dette i én fil
# LOADING MODIFIED PRIORS AND DATA
load(file = "original/mu_logit_p.RData")
load(file = "original/sigma_logit_p.RData")

# Loading data
load(file = "data/LDDdata.rda")

# READING THE MODEL CODE
source("R/nimble_models/nimbleCode_DOMM_lognormal_1.q")

# CONSTANTS USED FOR INITIAL INITUAL VALUES AND FOR PRIORS
colsumy = apply(LDDdata$data$y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
Nhat = apply(LDDdata$data$Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
lambdahat_surv = (Nhat/LDDdata$const$N_sites)/apply(LDDdata$const$area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)
lambdahat = tapply(lambdahat_surv, list(LDDdata$const$sam), mean)
N = round(LDDdata$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
N[is.na(N)] = 0
nrowY = nrow(LDDdata$data$Y)
ncolY = ncol(LDDdata$data$Y)

# FUNCTION FOR INITIAL VALUES
Inits = function(){
  sigma = runif(1, 0, 0.1)
  p1 = exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 = exp(log(p2hat)*runif(1, 0.9, 1.1))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    epsilon = matrix(rnorm(nrowY*ncolY, 0, sigma) , nrow = nrowY, ncol = ncolY),
    logit_p1 = rnorm(mu_logit_p, sigma_logit_p/5),
    logit_p2 = rnorm(mu_logit_p, sigma_logit_p/5),
    N = N,
    sigma = sigma,
    beta = runif(1, -0.5, 0.5)
  )
}

# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  nimbleCode_DOMM_lognormal_1,
  constants = list(lamblow = 0.1*lambdahat,  # 0.1 to 10 times point estimate
                   lambupp = 10*lambdahat,
                   N_surv = length(LDDdata$const$N_sites),
                   N_sites = LDDdata$const$N_sites,
                   mu_logit_p = mu_logit_p,
                   sigma_logit_p = sigma_logit_p,
                   sam = LDDdata$const$sam,
                   N_sam = length(unique(LDDdata$const$sam)),
                   area = LDDdata$const$area,
                   x = standardize(LDDdata$const$mean_field_dist)
  ),
  data = LDDdata$data,
  inits = Inits()
)

CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel) # Needs to be compiled for the last step
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("median_lambda", "mean_lambda" , "p1", "p2", "mu0", "sigma", "beta"), enableWAIC = TRUE)
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

posterior_lognormal <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=5000,#00,
  nburnin=1000,#00,
  nchain=3,
  thin=2,
  inits = Inits,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)

plot(posterior_lognormal$samples)
summary(posterior_lognormal$samples)
posterior_lognormal$WAIC

#save(posterior_lognormal, file = "Manuscript_models/posterior_lognormal.RData")
