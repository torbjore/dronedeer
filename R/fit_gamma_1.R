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

Inits <- function(){
  rate <- runif(1, 0.3, 1)
  p1 <- exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 <- exp(log(p2hat)*runif(1, 0.9, 1.1))
  list(
    #mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    exp_mu0 = lambdahat*runif(length(lambdahat), 0.9, 1.1),
    logit_p1 = rnorm(prior_parameters_for_p$mu_logit_p, prior_parameters_for_p$sigma_logit_p/5),
    logit_p2 = rnorm(prior_parameters_for_p$mu_logit_p, prior_parameters_for_p$sigma_logit_p/5),
    N = N,
    lambda = N + 0.01,
    invrate = 1/rate,
    beta = runif(1, -0.5, 0.5),
    New_Y = N
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

# dataNodes <- DoubleObsMultisiteModel$getNodeNames(dataOnly = TRUE)
# parentNodes <- DoubleObsMultisiteModel$getParents(dataNodes, stochOnly = TRUE, includeData = FALSE, upstream = TRUE)
# simNodes <- DoubleObsMultisiteModel$getDependencies(parentNodes, self = FALSE)

t1 <- Sys.time()
CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel) # Needs to be compiled for the last step
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, 
                                        monitors = c("p1", "p2", "mu0", "rate", "beta", "Disc_Y", "Disc_New_Y"), 
                                        #monitors = c("p1", "p2", "mu0", "rate", "beta"), 
                                        enableWAIC = TRUE)

# Setting up a block sampler (I've tried various combinations of blocking, chains often get stuck with any combination)
nn1 <- DoubleObsMultisiteModel$expandNodeNames(c('beta', 'rate')) #
DoubleObsMultisiteConf$removeSamplers(nn1)
DoubleObsMultisiteConf$addSampler(nn1, 'RW_block', control = list(adaptScaleOnly=FALSE))
# nn2 <- DoubleObsMultisiteModel$expandNodeNames(c('exp_mu0')) #
# DoubleObsMultisiteConf$removeSamplers(nn2)
# DoubleObsMultisiteConf$addSampler(nn2, 'RW_block', control = list(adaptScaleOnly=FALSE))

print(DoubleObsMultisiteConf)

DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

# DoubleObsMultisiteConf$getMonitors()
# DoubleObsMultisiteConf$getSamplerDefinition('beta')
# DoubleObsMultisiteConf$getSamplers('beta')
# DoubleObsMultisiteConf$printSamplers('exp_mu0')


t2 <- Sys.time()
cat("Compilation time:")
t2-t1

t3 <- Sys.time()
init.values <- list(Inits(), Inits(), Inits())
posterior_gamma_1 <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter = 10000, #115000,
  nburnin = 5000, #15000,
  nchain = 3,
  thin = 4,
  inits = init.values,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
t4 <- Sys.time()

cat("Run time:")
t4-t3

plot(posterior_gamma_1$samples) # 1 = svart, 2 = rød, 3 = grønn
summary(posterior_gamma_1$samples)
gelman.diag(posterior_gamma_1$samples)

lapply(posterior_gamma_1$samples, function(i) var(i[,"rate"]))
lapply(posterior_gamma_1$samples, function(i) mean(i[,"rate"]))
lapply(posterior_gamma_1$samples, function(i) mean(i[,"beta"]))
lapply(init.values, function(i) i[c("exp_mu0", "invrate", "beta")])

posterior_gamma_1$WAIC

# nimbleList object of type waicNimbleList
# Field "WAIC":
# ORIGINAL  [1] 117.1246
#   Field "WAIC":
# [1] 117.1465
#   


#save(posterior_gamma_1, file = "data/posterior_samples/gamma_1.RData")

# # Posterior predictive checks
# 
# samp <- as.matrix(posterior_gamma_1$samples)
# plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"])
# abline(0,1, col="red")
# 
# mean(samp[,"Disc_New_Y"] > samp[,"Disc_Y"])

