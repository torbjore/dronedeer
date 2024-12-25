# R-packages used: 
library(nimble)
library(coda)

# Making standardize function (NB! Do not use scale() as this function scales each column of a matrix)
standardize <- function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

# Loading prior-parameters for p
load(file = "data/prior_parameters_for_p.rda")

# FOR WIDE p
prior_parameters_for_p$mu_logit_p <- 0
prior_parameters_for_p$sigma_logit_p <- 2

# Loading data
load(file = "data/LDDdata.rda")

# READING THE MODEL CODE
source("R/nimble_models/nimbleCode_DOMM_gamma_2_with_N.q")

# CONSTANTS USED FOR INITIAL INITUAL VALUES AND FOR PRIORS
colsumy <- apply(LDDdata$data$y, 3, sum, na.rm=TRUE)
p1hat <- colsumy[3]/(colsumy[2]+colsumy[3])
p2hat <- colsumy[3]/(colsumy[1]+colsumy[3])
Nhat <- apply(LDDdata$data$Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
# lambdahat_surv <- (Nhat/LDDdata$const$N_sites)/apply(LDDdata$const$area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)
# lambdahat <- tapply(lambdahat_surv, list(LDDdata$const$sam), mean)
lambdahat <- LDDdata$const$lambdahat

N <- round(LDDdata$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
N[is.na(N)] <- 0
nrowY <- nrow(LDDdata$data$Y)
ncolY <- ncol(LDDdata$data$Y)

Inits <- function(){
  sigma_p <- runif(1, 0.1, 0.4)
  p1 <- exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 <- exp(log(p2hat)*runif(1, 0.9, 1.1))
  mu_p1 <- log(p1/(1-p1))
  mu_p2 <- log(p2/(1-p2))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    mu_p1 = mu_p1,
    mu_p2 = mu_p2,
    logit_p1 = matrix(rnorm(nrowY*ncolY, mu_p1, sigma_p/2) , nrow = nrowY, ncol = ncolY),
    logit_p2 = matrix(rnorm(nrowY*ncolY, mu_p2, sigma_p/2) , nrow = nrowY, ncol = ncolY),
    N = N,
    lambda = N + 0.01,

    beta = runif(2, -0.5, 0.5),
    sigma = runif(1, 0.5, 1),
    New_Y = N, # Warning message if not included
    New_y = LDDdata$data$y
  )
}

# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  nimbleCode_DOMM_gamma_2,
  constants = list(#lamblow = 0.1*lambdahat,  # 0.1 to 10 times point estimate
                   #lambupp = 10*lambdahat,
                   #mu0hat = log(lambdahat/exp(0.5*1.47^2)),
                   N_surv = length(LDDdata$const$N_sites),
                   N_sites = LDDdata$const$N_sites,
                   prior_mu_logit_p = prior_parameters_for_p$mu_logit_p,
                   prior_sigma_logit_p = prior_parameters_for_p$sigma_logit_p,
                   sam = LDDdata$const$sam,
                   N_sam = length(unique(LDDdata$const$sam)),
                   S = sapply(1:8, function(i) range(which(LDDdata$const$sam==i))),
                   area = LDDdata$const$area,
                   x = standardize(LDDdata$const$mean_field_dist)
  ),
  data = LDDdata$data,
  inits = Inits()
)

t1 <- Sys.time()
CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel) # Needs to be compiled for the last step
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, 
                                        monitors = c("Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "mu_p1", "mu_p2", "mu0", "sigma", "sigma_p", "beta", "Dens_tot"),
                                        enableWAIC = TRUE)

# Setting up a block sampler (I've tried various combinations of blocking, chains often get stuck with any combination)
# nn1 <- DoubleObsMultisiteModel$expandNodeNames(c('beta', 'rate')) #
# DoubleObsMultisiteConf$removeSamplers(nn1)
# DoubleObsMultisiteConf$addSampler(nn1, 'RW_block', control = list(adaptScaleOnly=FALSE))
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
settings <- list(
  niter = 200000,
  nburnin = 20000,
  nchain = 3,
  thin = 6
  # # FOR LONG RUN
  # niter = 1100000,
  # nburnin = 20000,
  # nchain = 3,
  # thin = 18
)
out <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter = settings$niter,
  nburnin = settings$nburnin,
  nchain = settings$nchain,
  thin = settings$thin,
  inits = init.values,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
t4 <- Sys.time()
cat("Run time:")
t4-t3

# Saving workspace
#save(settings, out, file = "data/posterior_samples/gamma_2_LONGrun_with_N.RData")
#save(settings, out, file = "data/posterior_samples/gamma_2_with_N_wide_p.RData")

load("data/posterior_samples/gamma_2_with_N_wide_p.RData")
out_w <- out
load("data/posterior_samples/gamma_2_with_N.RData")

#plot(out$samples) # 1 = black, 2 = red, 3 = green
summary(out$samples)
summary(out_w$samples)

# gelman.diag(out$samples)

# out$WAIC
# 
# # # Posterior predictive checks
# # 
# samp <- as.matrix(out$samples)
# plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"])
# abline(0, 1, col="red")
# mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])
# # OK!
# 
# # Wrt y
# plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
# abline(0, 1, col="red")
# mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])

