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
source("R/nimble_models/nimbleCode_DOMM_gamma_1.q")

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

# MERK: Ser at vi ikke har noen starting values for lambda

Inits = function(){
  rate = runif(1, 0.5, 5)
  p1 = exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 = exp(log(p2hat)*runif(1, 0.9, 1.1))
  list(
#    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    exp_mu0 = lambdahat*runif(length(lambdahat), 0.9, 1.1),
    logit_p1 = rnorm(prior_parameters_for_p$mu_logit_p, prior_parameters_for_p$sigma_logit_p/5),
    logit_p2 = rnorm(prior_parameters_for_p$mu_logit_p, prior_parameters_for_p$sigma_logit_p/5),
    N = N,
    rate = rate,
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
#DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("mean_lambda" , "p1", "p2", "mu0", "rate", "beta"), enableWAIC = TRUE)
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("p1", "p2", "exp_mu0", "rate", "beta"), enableWAIC = TRUE)
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

# Trying to set up a block sampler for 'rate', 'mu0[1:8]' nodes
nn <- DoubleObsMultisiteModel$expandNodeNames(c('rate', 'beta', 'mu0'))
DoubleObsMultisiteConf$removeSamplers(nn)
DoubleObsMultisiteConf$addSampler(nn, 'RW_block', control = list(adaptScaleOnly=FALSE))
#DoubleObsMultisiteConf$addSampler(nn, 'AF_slice', control = list(adaptScaleOnly=FALSE))

t2 <- Sys.time()
cat("Compilation time:")
t2-t1

t3 <- Sys.time()
posterior_gamma_1 <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter= 6000, #0,
  nburnin= 5000,
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

lapply(posterior_gamma_1$samples, function(i) apply(i, 2, mean)[c("rate", "beta")])
# # Run 1:
# $chain1
# rate     beta 
# 0.162997 0.322608 
# 
# $chain2
# rate      beta 
# 0.1563404 0.6641122 
# 
# $chain3
# rate       beta 
# 0.05897441 0.42049683


posterior_gamma_1$WAIC

# nimbleList object of type waicNimbleList
# Field "WAIC":
# ORIGINAL  [1] 117.1246
#   Field "WAIC":
# [1] 117.1465
#   

gelman.diag(posterior_gamma_1$samples)

#save(posterior_gamma_1, file = "data/posterior_samples/gamma_1.RData")