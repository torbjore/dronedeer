# R-packages used: 
library(nimble)
library(coda)

# Making standardize function (NB! Do not use scale() as this function scales each column of a matrix)
standardize <- function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

# Loading prior-parameters for p
load(file = "data/prior_parameters_for_p.rda")

# Loading data
load(file = "data/LDDdata.rda")

# READING THE MODEL CODE
source("R/nimble_models/nimbleCode_DOMM_gamma_1_original.q")

# CONSTANTS USED FOR INITIAL INITUAL VALUES AND FOR PRIORS
colsumy <- apply(LDDdata$data$y, 3, sum, na.rm=TRUE)
p1hat <- colsumy[3]/(colsumy[2]+colsumy[3])
p2hat <- colsumy[3]/(colsumy[1]+colsumy[3])
Nhat <- apply(LDDdata$data$Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
lambdahat_surv <- (Nhat/LDDdata$const$N_sites)/apply(LDDdata$const$area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)
lambdahat <- tapply(lambdahat_surv, list(LDDdata$const$sam), mean)
N <- round(LDDdata$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
N[is.na(N)] <- 0
nrowY <- nrow(LDDdata$data$Y)
ncolY <- ncol(LDDdata$data$Y)

Inits = function(){
  sigma = runif(1, 1, 2)
  p1 = exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 = exp(log(p2hat)*runif(1, 0.9, 1.1))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    logit_p1 = rnorm(prior_parameters_for_p$mu_logit_p, prior_parameters_for_p$sigma_logit_p/5),
    logit_p2 = rnorm(prior_parameters_for_p$mu_logit_p, prior_parameters_for_p$sigma_logit_p/5),
    N = N,
    sigma = sigma,
    beta = runif(1, -0.5, 0.5)
  )
}

# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  nimbleCode_DOMM_gamma_1,
  constants = list(lamblow = 0.1*lambdahat,  # 0.1 to 10 times point estimate
                   lambupp = 10*lambdahat,
                   N_surv = length(LDDdata$const$N_sites),
                   N_sites = LDDdata$const$N_sites,
                   mu_logit_p = prior_parameters_for_p$mu_logit_p,
                   sigma_logit_p = prior_parameters_for_p$sigma_logit_p,
                   sam = LDDdata$const$sam,
                   N_sam = length(unique(LDDdata$const$sam)),
                   area = LDDdata$const$area,
                   x = standardize(LDDdata$const$mean_field_dist)
  ),
  data = LDDdata$data,
  inits = Inits()
)

t1 <- Sys.time()
CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel) # Needs to be compiled for the last step
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("mean_lambda" , "p1", "p2", "mu0", "sigma", "beta"), enableWAIC = TRUE)
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)
t2 <- Sys.time()

cat("Compilation time:")
t2-t1

t3 <- Sys.time()
posterior_gamma_1 <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=5000, #0,
  nburnin=1000,#0,
  nchain=3,
  thin=4,
  inits = Inits,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
t4 <- Sys.time()

cat("Run time:")
t4-t3

plot(posterior_gamma_1$samples)
summary(posterior_gamma_1$samples)
posterior_gamma_1$WAIC

# nimbleList object of type waicNimbleList
# Field "WAIC":
#   [1] 117.1246

gelman.diag(posterior_gamma_1$samples)

#save(posterior_gamma_1, file = "data/posterior_samples/gamma_1.RData")