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
source("R/nimble_models/nimbleCode_DOMM_exponential_1.q")

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

# FUNCTION FOR INITIAL VALUES
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
    sigma_p = sigma_p,
    beta = runif(1, -0.5, 0.5),
    New_Y = N, # Warning message if not included
    New_y = LDDdata$data$y
  )
}

# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  nimbleCode_DOMM_exponential_1,
  constants = list(lamblow = 0.1*lambdahat,  # 0.1 to 10 times point estimate
                   lambupp = 10*lambdahat,
                   N_surv = length(LDDdata$const$N_sites),
                   N_sites = LDDdata$const$N_sites,
                   prior_mu_logit_p = prior_parameters_for_p$mu_logit_p,
                   prior_sigma_logit_p = prior_parameters_for_p$sigma_logit_p,
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
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, 
                                        monitors = c("Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "mu_p1", "mu_p2", "mu0", "sigma_p", "beta"), enableWAIC = TRUE)


DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)
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
save(settings, out, file = "data/posterior_samples/exponential_1_run1.RData")

#plot(out$samples) # 1 = black, 2 = red, 3 = green
# summary(out$samples)
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
